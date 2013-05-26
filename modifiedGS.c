#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<omp.h>
#include<math.h>
#include<assert.h>
#include<pthread.h>
#include<limits.h>
#include<semaphore.h>
#include<float.h>
#include<sys/time.h>

#define OMP_PTHREADS 1

#define DEBUG 0
int should_die = 0;
sem_t thread_sema;
int num_threads = 8;
pthread_t match_thread[24];
int thread_ids[24];
int thread_data[24]; /* Int Array for thread data. USAGE: memset in your function and access with thread_data[tid]. 
		      * Create new array for a different data type
		      */
pthread_attr_t attr;
pthread_mutex_t q_mutex;
int current_color = 0;

/*
 * Structure to hold an edge.
 * Single long long integer is used to hold two integers
 * Given an edge (a , b), assuming both integers:
 * Last 32 bits hold 'b' and the prev 32 bits hold 'a'
 */

typedef struct {
    long long int x;
    pthread_mutex_t e_mutex;
    int weight; /* FIXME: should be float? */
    int valid;
} edge;

typedef struct {
    int id;/* index of the neighbor */
    int weight;/* weight of the edge to the  neighbor */
} neighbor;

/*
 * Structure to hold vertex info.
 * degree of vertex
 * index of vertex is implicilty stored in the accessing array iterator
 */

typedef struct {
    int id;/* integer index */
    int deg;/* degree of vertex */
    int inC; 
    int inL;
    int inS;/* is in the current MIS */
    int inW;/* to be considered for next MIS */
    int color; 
    int match; /* if matched, its match; else -1 */
    neighbor *neighbors; /* adjacency list */
    int neighbor_count; 
    pthread_mutex_t v_mutex;
    int weight;/* FIXME: should be float? */
    int predecessor_id[2]; /* This vertex's association with atmost 2 vertices in previous graph */
    int successor_id;   /* This vertex's association with the vertex in next graph */
} vertex;

typedef struct graph {
    pthread_mutex_t g_mutex;
    int num_v; /* number of vertices */
    int num_e; /* number of edges */
    vertex *vl; /* vertex list */
    edge *edl;  /* edge list */
    int uncolored; /* number of vertices not colored */
    int unmatched; /* number of vertices not matched */
    int num_colors;
    int max_degree;
    struct graph *next;/* pointer to next graph */
} graph;

typedef void*(*funct)(void*);
/* element of the work queue */
typedef struct work_elem {
  funct func;
  int bucket;
  struct work_elem *next;
} work;

typedef struct work_queue {
  int size;
  work *head;
  work *tail;
  pthread_mutex_t wq_mutex; // for mutual exclusion while adding/removing element
  sem_t wq_sema;
} work_queue;

typedef struct queue_elem {
    int v; /* index of the vertex */
    struct queue_elem *next;
} elem;

typedef struct queue {
    int size;
    elem *head;
} queue;

graph *gr_root;
graph *cur_gr;/* Points to the current graph being processed. Used by threads for sharing */
graph *next_level_gr;/* Points to the graph next to current graph being processed. Only used in build_adjacnecy_list_next_graph() */
queue *q;
work_queue* wqueue[24];
pthread_t thread[24];

static queue* queue_init() {
    queue *q = malloc(sizeof(queue));
    if (!q) {
	printf("%s: malloc failed \n", __func__);
	exit(-1);
    }
    q->size = 0;
    q->head = NULL;
    pthread_mutex_init(&q_mutex, NULL);
    return q;
}

work_queue* wq_init() {
  work_queue *wq = malloc(sizeof(work_queue));
  if (!wq) {
    printf("%s: malloc failed \n", __func__);
    exit(-1);
  }
  wq->size = 0;
  wq->head = NULL;
  wq->tail = NULL;
  pthread_mutex_init(&wq->wq_mutex, NULL);
  if (sem_init(&wq->wq_sema, 0, 0)) {
    printf("sem_init failed\n");
    exit(0);
  }
  return wq;
}

static void queue_free(queue *q) {
    int i;
    assert(q);
    elem *cur = q->head;
    elem *temp;
    while(cur) {
	temp = cur->next;
	free (cur);
	cur = temp;
    }
    free(q);
    pthread_mutex_destroy(&q_mutex);
}

static void queue_add(queue *q, int v) {
    elem *e = malloc(sizeof(elem));
    if (!e) {
	printf("%s: malloc failed\n", __func__);
	exit(-1);
    }
    e->v = v;
    /* add element into the front of the queue */
    pthread_mutex_lock(&q_mutex);
    e->next = q->head;
    q->head = e;
    q->size++;
    pthread_mutex_unlock(&q_mutex);
}

void wq_add(work_queue *wq, funct f) {
  work *w = calloc(1, sizeof(work));
  if (!w) {
    printf("%s: malloc failed\n", __func__);
    exit(-1);
  }
  w->func = f;
  w->next = NULL;
  pthread_mutex_lock(&wq->wq_mutex);

  /* add element at the end of queue */
  if (wq->head == NULL) {
    wq->head = w;
    wq->tail = w;
  } else {
    wq->tail->next = w;
    wq->tail = w;
  }
  wq->size++;
  pthread_mutex_unlock(&wq->wq_mutex);
#if DEBUG > 1
  printf("adding an element to wq: size=%d\n", wq->size);
#endif
  sem_post(&wq->wq_sema);
}

work* wq_remove(work_queue *wq) {

  work *ret;
  /* remove element from the front of the queue */
  pthread_mutex_lock(&wq->wq_mutex);
#if DEBUG > 1
  printf("%s: size=%d, head = %p, tail = %p\n", __func__, wq->size, (void *)wq->head, (void *)wq->tail);
#endif
  if (wq->head == NULL || wq->size == 0)
    return NULL;

  ret = wq->head;
  wq->head = ret->next;
  wq->size--;
  if (wq-> size == 1) {
    assert(wq->head == wq->tail);
  }
  if(wq->head == NULL) {
    wq->tail = NULL;
  }
  pthread_mutex_unlock(&wq->wq_mutex);
  return ret;
}

void *thread_func(void *id) {

  int tid = *(int *)(id);
  work *w;
  int value;
  while (1) {

#if DEBUG > 1
    printf("Thread[%d]: going to sleep...\n", tid);
#endif
    //sem_getvalue(&wqueue[tid]->wq_sema, &value);
    //printf("The value of the semaphore is %d\n", value);
    sem_wait(&wqueue[tid]->wq_sema);

    if (should_die) {
      printf("Thread[%d] exiting...\n", tid);
      pthread_exit(NULL);
    }

    /* 1. wait till an element is in work_queue
           2. pick a work from work_queue
           3. Execute it
    */
#if DEBUG > 1
    printf("Thread[%d]: woke up\n", tid);
#endif
    w = wq_remove(wqueue[tid]);
#if DEBUG > 1
    printf("Thread[%d]: died in wq_remove\n", tid);
#endif
    assert(w);
    assert(w->func);
#if DEBUG > 1
    printf("Thread[%d]: died in assert\n", tid);
#endif
    w->func(id);
    free(w);
  }
}

void safe_decrement(int *x) {
    graph *gr = cur_gr;
    pthread_mutex_lock(&gr->g_mutex);
    *x = *x - 1;
    pthread_mutex_unlock(&gr->g_mutex);
}

void safe_increment(int *x) {
    graph *gr = cur_gr;
    pthread_mutex_lock(&gr->g_mutex);
    *x = *x + 1;
    pthread_mutex_unlock(&gr->g_mutex);
}

void setup_threads(int count)
{
  int i;
  int ret;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  /* initialize thread ids */
  for (i = 0; i < count; i++) {
    thread_ids[i] = i;
  }

  for (i = 0; i < count; i++) {
    ret = pthread_create(&thread[i], &attr, thread_func, (void *) &thread_ids[i]);
    if (ret) {
      printf("Error creating thread %d : %d \n", i, ret);
      exit(-1);
    } else {
#if DEBUG > 1
      printf("created thread %d\n", i);
#endif
    }
  }

  /* initialize semaphores */
  sem_init(&thread_sema, 0, 0);
}


/* create requested number of joinable threads */
static void create_threads(int count, void* (*work) (void *))
{
    int i;
    int ret;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
    for (i = 0; i < count; i++) {
	ret = pthread_create(&match_thread[i], &attr, work, (void *) &thread_ids[i]);
	if (ret) {
	    printf("Error creating thread %d : %d \n", i, ret);
	    exit(-1);
	} else {
	    #if DEBUG > 1
	    printf("created thread %d\n", i);
	    #endif
	}
    }
}


/* wait for all the threads to complete */
static void join_threads(int count) 
{
    int i, ret;
    void *status;

    for (i = 0; i < count; i++) {
	ret = pthread_join(match_thread[i], &status);
	if (ret) {
	    printf("Error joining thread %d : %d\n", i, ret);
	} else {
	    #if DEBUG > 1
	    printf("joining thread %d\n", i);
	    #endif
	}
    }
}


static void printinSVerts(vertex *vl, int numvertices)
{
    int i = 0;

    for ( i = 0; i < numvertices; i++) {
	if (vl[i].inS) {
	    printf("Vertex = %d Degree = %d\n", i, vl[i].deg);
	}
    }
}

int isX_greatestNeighborOf_Y_inS(vertex *vl, int ix, int iy)
{
    int i = 0, id = 0;
    for (i = 0; i < vl[iy].neighbor_count; i++) {
	id = vl[iy].neighbors[i].id;
	if (vl[id].inS) {
	    if (vl[ix].deg < vl[id].deg) {
		#if DEBUG > 1
		printf("vertex = %d not greatest neighbor of vertex = %d in S\n", ix, iy);
		#endif
		return 0;
	    }
	    else if (vl[ix].deg == vl[id].deg) {
		if (vl[ix].weight < vl[id].weight) {
		    #if DEBUG > 1
		    printf("vertex = %d not greatest neighbor of vertex = %d in S\n", ix, iy);
		    #endif
		    return 0;
		}
	    }
	}
    }

#if DEBUG > 1
    printf("vertex = %d is the greatest neighbor of vertex = %d in S\n", ix, iy);
#endif
    return 1;
}

/* Checks first for greatest degree and if not then for greatest weight among neighbors */
int greatestAmongNeighbors(vertex *vl, int id)
{
    int i = 0, greatestDeg = 1, greatestWeight = 1, idx = 0;
    /*
     * 1. If my degree greater than any of my neighbors's degree, return true.
     * 2. If my degree less than or equal to my neighbors's degree, check if
     *    my weight greater than any of my neighbor's weight. If yes return true, else false.
     */

    for (i = 0; i < vl[id].neighbor_count; i++) {
	idx = vl[id].neighbors[i].id;
	if (vl[idx].inL) {
	    if (vl[id].deg < vl[idx].deg) {
		return 0;
	    }
	    else if (vl[id].deg == vl[idx].deg) {
		if (vl[id].weight < vl[idx].weight) {
		    return 0;
		}
	    }
	}
    }

    return 1;
}

/*
 * FIXME: Walk the graph linked list starting from graph_root and 
 * free all graphs
 */
static void free_graph(graph *gr) {
    int i;
    assert(gr);
    for (i = 0; i < gr->num_v; i++) {
	free(gr->vl[i].neighbors);
    }
    if (gr->vl) free(gr->vl);
    if (gr->edl) free(gr->edl);
    pthread_mutex_destroy(&gr->g_mutex);
    free(gr);
}

static graph* create_graph() {
    graph *gr = malloc(sizeof(graph));
    if (!gr) {
	printf("%s: malloc failed\n", __func__);
	exit(-1);
    }

    gr->vl = NULL;
    gr->edl = NULL;
    gr->max_degree = 0;
    gr->num_colors = 0;
    gr->uncolored = 0;
    gr->unmatched = 0; 
    pthread_mutex_init(&gr->g_mutex, NULL);
    gr->next = NULL;

    return gr;
}

static void print_full_graph(graph *gr, int match)
{
    int i;
    for (i = 0; i < gr->num_v; i++) {
	int j;

	if (match) {
	    /* print match in the parantheses */
	    printf("%d(match:%d)", i, gr->vl[i].match);
	} else {
	    /* print color in the parantheses */
	    printf("%d(color:%d)", i, gr->vl[i].color);
	}
	printf(" (sid:%d) ", gr->vl[i].successor_id);
	printf("(pid[0]:%d pid[1]:%d)\n", gr->vl[i].predecessor_id[0], gr->vl[i].predecessor_id[1]);

	for(j = 0; j < gr->vl[i].neighbor_count; j++) {
	    printf("->");
	    printf("%d(%d)", gr->vl[i].neighbors[j].id, gr->vl[i].neighbors[j].weight);
	}
	printf("\n");
    }
}

static void populate_graph(graph *gr, char *file) 
{
    
    FILE *fp;
    int i = 0, x, y, w;
    int j = 0;
    long long e;
    
    fp = fopen(file, "r");
    fscanf(fp, "%d %d", &gr->num_v, &gr->num_e);
    gr->uncolored = gr->num_v;
    gr->unmatched = gr->num_v; 
    
    /* I assumed that there are more vertices than threads.
     * If this assert hits, change the graph or num_threads.
     * The issue will be with sharing across threads.
     */
    assert(gr->num_v > num_threads);
    
    gr->edl = (edge *)calloc(gr->num_e, sizeof(edge));
    if (!gr->edl) {
	printf("Failed to allocate memory for graph edges. Bailing out.\n");
	exit(0);
    }
    
    gr->vl = (vertex *)calloc(gr->num_v, sizeof(vertex));
    if (!gr->vl) {
	printf("Failed to allocate memory for graph vertices. Bailing out.\n");
	free_graph(gr);
	exit(0);
    }
    
    for (i = 0; i < gr->num_e; i++) {
	fscanf(fp, "%d %d %d", &x, &y, &w);
	e = 0;
	x -= 1;
	y -= 1;
	e = (long long int)x;
	e = e << 32;
	e = e | (long long int)y;
	gr->edl[i].x = e;
	gr->edl[i].valid = 1;
	gr->edl[i].weight = w;
	#if DEBUG > 1
	printf("Added edge %d,%d\n", x, y);
	#endif
	
	/* for vertex 'x' */
	gr->vl[x].id = x;
	gr->vl[x].deg++;
	gr->vl[x].color = 0; /* colors start from 1 */
	gr->vl[x].match = -1; /* vertex is not yet matched */
	gr->vl[x].inC = 1;
	gr->vl[x].inS = 0;
	gr->vl[x].inL = 0;
	gr->vl[x].predecessor_id[0] = -1;
	gr->vl[x].predecessor_id[1] = -1;
	gr->vl[x].successor_id = -1;
	if (gr->max_degree < gr->vl[x].deg) { gr->max_degree = gr->vl[x].deg; }
	
	/* for vertex 'y' */
	gr->vl[y].id = y;
	gr->vl[y].deg++;
	gr->vl[y].color = 0;
	gr->vl[y].match = -1;
	gr->vl[y].inC = 1;
	gr->vl[y].inS = 0;
	gr->vl[y].inL = 0;
	gr->vl[y].predecessor_id[0] = -1;
	gr->vl[y].predecessor_id[1] = -1;
	gr->vl[y].successor_id = -1;
	if (gr->max_degree < gr->vl[y].deg) { gr->max_degree = gr->vl[y].deg; }
    }
    
    /* Loop through the edges to create adjacency list for each vertex */
    for (i = 0; i < gr->num_e; i++) {
	x = gr->edl[i].x >> 32;
	y = gr->edl[i].x & 0xFFFFFFFF;
	if (gr->vl[x].neighbors == NULL) {
	    gr->vl[x].neighbors = (neighbor *)calloc(gr->max_degree, sizeof(neighbor));
	}
	if (gr->vl[y].neighbors == NULL) {
	    gr->vl[y].neighbors = (neighbor *)calloc(gr->max_degree, sizeof(neighbor));
	}
	gr->vl[x].neighbors[gr->vl[x].neighbor_count].id = y;
	gr->vl[x].neighbors[gr->vl[x].neighbor_count].weight = gr->edl[i].weight;
	gr->vl[x].neighbor_count++;
	
	gr->vl[y].neighbors[gr->vl[y].neighbor_count].id = x;
	gr->vl[y].neighbors[gr->vl[y].neighbor_count].weight = gr->edl[i].weight;
	gr->vl[y].neighbor_count++;
    }
}

/* get the smallest color not taken by any 
 *  of the neighbors of the given vertex id
 */
static int get_best_color(int id)
{
    graph *gr = cur_gr;
    int i, j, idx;
    int min = 1, max = gr->num_colors;
    int ret = 1;
    int is_seen = 0;
    
    assert(gr->vl[id].color < 1); /* shouldn't be colored earlier */
    for (i = min; i <= max; i++) {
	is_seen = 0;
	for (j = 0; j < gr->vl[id].neighbor_count; j++) {
	    idx = gr->vl[id].neighbors[j].id;
	    if (gr->vl[idx].color == i) {
		is_seen = 1;
	    }
	}
	if (is_seen == 0) 
	    return i;
    }
    return i;
}

/* find the neighbor for vertex 'src'
 * that is not yet matched and has the
 * highest weight.
 * FIXME: for now weight is just the index
 * New Policy: If my index is greatest of all my neighbors, go ahead else quit.
 * Also find the greatest indexed one of all my neighbors. Then pick that neighbor and see if I am its greatest indexed neighbor.
 * If true pair myself with that neighbor.
 */
int get_best_unmatched_neighbor(graph *gr, int src) {
    
    int i = 0, ret = -1, max_w = -1, idx = 0;
    for (i = 0; i < gr->vl[src].neighbor_count; i++) {
	idx = gr->vl[src].neighbors[i].id;
#if 0
	if (gr->vl[idx].match == -1 && gr->vl[idx].weight > max_w) {
	    max_w = gr->vl[idx].weight;
	    ret = idx;
	}
#endif
	if (gr->vl[idx].match == -1 && idx > max_w) {
	    max_w = idx;
	}
    }
    return max_w;
#if 0
    int greatest_neighbor = -1;
    int idx, i;
    for (i = 0; i < gr->vl[src].neighbor_count; i++) {
	idx = gr->vl[src].neighbors[i].id;
	if (idx > greatest_neighbor) greatest_neighbor = idx;
    }
    if (src < greatest_neighbor) {
	return -1;/* let the greatest neighbor handle matching */
    }

    /* Check if I am the greatest neighbor of my greatest neighbor */
    for (i = 0; i < gr->vl[greatest_neighbor].neighbor_count; i++) {
	idx = gr->vl[greatest_neighbor].neighbors[i].id;
	if (idx > src) return -1;
    }

    return greatest_neighbor;
#endif
}

void* do_matching(void *tid) {

    graph *gr = cur_gr;
    int share = (gr->num_v)/num_threads;
    int start = *(int *)(tid) * share;
    int end = start + share - 1;
    int i;
    int src, target; //(u, v) in algo

    /* last threads gets all remaining */
    if (*(int *)tid == num_threads - 1) {
	end = gr->num_v - 1;
    }

#if DEBUG > 1
    printf("%s: I am thread[%d] working on [%d to %d]\n", __func__, *(int*)tid, start, end);
#endif

    /* no locking needed because each thread
     *  get mutually exclusive share of the vertices
     */
    for (i = start; i <= end; i++) {
	src = i;
	/* unmatched vertices of this color */
	if (gr->vl[src].color == current_color && 
	    gr->vl[src].match == -1) {
	    /* FIXME: ensures ties are broken and only one thread changes the match */
	    target = get_best_unmatched_neighbor(gr, src);
	    if ( target > -1) {
		gr->vl[src].match = target;
		gr->vl[target].match = src;
		safe_decrement(&gr->unmatched);
		safe_decrement(&gr->unmatched); 
		queue_add(q, src);
	    }
	}
    }

    pthread_attr_destroy(&attr);
    pthread_exit(NULL);
}

void maximal_matching(graph *gr)
{
    int colors = gr->num_colors;
    int i;
    elem *cur;
    int src_match, dst_match, dst, j;

    for (i = 1; i <= gr->num_colors; i++) {
	current_color = i;
	q = queue_init();

#if 0
	#pragma omp for
	for (j = 0; j < gr->num_v; j++) {
	    src = j;
	    /* unmatched vertices of this color */
	    if (gr->vl[src].color == i && 
		gr->vl[src].match == -1) {
		/* FIXME: ensures ties are broken and only one thread changes the match */
		target = get_best_unmatched_neighbor(gr, src);
		gr->vl[src].match = target;
		gr->vl[target].match = src;
		queue_add(q, src);
	    }
	}
#endif
	create_threads(num_threads, do_matching);
	join_threads(num_threads);
	
	cur = q->head;

	while(cur) {
	    src_match = gr->vl[cur->v].match;
	    dst = src_match;
	    dst_match = gr->vl[dst].match;
	    if (src_match > -1 && dst_match != cur->v) {
		gr->vl[cur->v].match = -1;
		safe_increment(&gr->unmatched);
		safe_increment(&gr->unmatched);
	    }
	    cur = cur->next;
	}

	if (((gr->num_v - gr->unmatched) % 2) != 0) {
	    printf("======= Error in matching Total vertices %d, unmatched %d====== \n", gr->num_v, gr->unmatched);
	    #if 0
	    print_full_graph(cur_gr, 1);
	    #endif
	    for (i = 0; i < gr->num_v; i++) {
		if (gr->vl[i].match == -1) {
		    printf("Unmatched %d(-1)\n", i);
		} else {
		    if (i != gr->vl[gr->vl[i].match].match) {
			printf("Incorrect match %d matched to (%d) while %d matched to (%d)", i, gr->vl[i].match, gr->vl[i].match, gr->vl[gr->vl[i].match].match);
		    } else {
			printf("Correct match %d(%d)\n", i, gr->vl[i].match);
		    }
		}
	    }
	}
	assert(((gr->num_v - gr->unmatched) % 2) == 0);
	queue_free(q);
    }
}

void * mis_transition_from_L_to_S(void *thread_id)
{
    graph *gr = cur_gr;
    int share = (gr->num_v)/num_threads;
    int tid = *(int *)(thread_id);
    int start = tid * share;
    int end = start + share - 1;
    int i;

    /* last threads gets all remaining */
    if (tid == num_threads - 1) {
	end = gr->num_v - 1;
    }

#if DEBUG > 1
    printf("%s: I am thread[%d] working on [%d to %d]\n", __func__, tid, start, end);
#endif

    /* no locking needed because each thread
     *  get mutually exclusive share of the vertices
     */
    for (i = start; i <= end; i++) {
	if(gr->vl[i].inL) {
	    if(greatestAmongNeighbors(gr->vl, i)) {
		/* Put yourself in S. Remove yourself from C */
		if(gr->vl[i].inS) {
		    continue;
		}
		gr->vl[i].inS = 1;
		gr->vl[i].inC = 0;
		#if DEBUG > 1
		printf("Putting vertex = %d  weight = %d in S\n", i, gr->vl[i].weight);
		#endif
		thread_data[tid]++;
	    }
	}
    }

    pthread_attr_destroy(&attr);
    pthread_exit(NULL);
}

void * mis_transition_remove_neighbors_of_S_from_C(void *thread_id)
{
    graph *gr = cur_gr;
    int share = (gr->num_v)/num_threads;
    int tid = *(int *)(thread_id);
    int start = tid * share;
    int end = start + share - 1;
    int i, j;

    /* last threads gets all remaining */
    if (tid == num_threads - 1) {
	end = gr->num_v - 1;
    }

#if DEBUG > 1
    printf("%s: I am thread[%d] working on [%d to %d]\n", __func__, tid, start, end);
#endif
    /* no locking needed because each thread
     *  get mutually exclusive share of the vertices
     */
    for (i = start; i <= end; i++) {
	if(gr->vl[i].inS) {
	    gr->vl[i].inL = 0;
	    /* Walk through my neighbor list. Delete each neighbor if I am its greatest neighbor */
	    for (j = 0; j < gr->vl[i].neighbor_count; j++) {
		int idx = gr->vl[i].neighbors[j].id;
		#if DEBUG > 1
		printf("vertex = %d neighbor count = %d  Neighbor vertex = %d\n", i, gr->vl[i].neighbor_count, idx);
		#endif
		if (gr->vl[idx].inC) {
		    if(isX_greatestNeighborOf_Y_inS(gr->vl, i, idx)) {
			gr->vl[idx].inC = 0;
			gr->vl[idx].inL = 0;
			thread_data[tid]++;
			#if DEBUG > 1
			printf("vertex = %d  removed vertex = %d from C\n", i, idx);
			#endif
		    }
		}
	    }
	}
    }

    pthread_attr_destroy(&attr);
    pthread_exit(NULL);
}

static void omp_mis_luby(graph *gr)
{
    int i = 0, x, y, nvertsinC = 0, nedges = 0;
    int j = 0, itr = -1;
    long long maxweight, e;
    float sum_of_sel_probability = 0.0;


    /* Run luby's only for uncolored vertices */
    for (i = 0; i < gr->num_v; i++) {
	if(gr->vl[i].color == 0) { /* not yet colored */
	    gr->vl[i].inC = 1; /* include them in the MIS */
	    nvertsinC++;
	}
    }
    
    maxweight = pow(nvertsinC, 4);
    
    while (nvertsinC) {
	int sum_removed = 0, num_verts_to_choose = 0;
	sum_of_sel_probability = 0.0;
	itr++;
	
	// printf("Iteration = %d vertices in C = %d\n", itr, nvertsinC);
	
	/* Clear thread_data before usage */
	memset(thread_data, 0, sizeof(thread_data));
	
	for (i = 0; i < gr->num_v; i++) {
	    if(gr->vl[i].inC) {
		gr->vl[i].weight = rand() % maxweight;
		sum_of_sel_probability += 1.0/(2 * (gr->vl[i].deg));
	    }
	}
	
	/* Choose a random number k between 2 to numvertices/2 
	 * Select 1 vertex based on 1/2d probability
	 * Select a random weight between 0 to sum_of_sel_probability
	 * Find the weight that causes current sum to cross the random number
	 * Repeat k times
	 */
	num_verts_to_choose = nvertsinC; 
	if(nvertsinC > num_threads) {
	    num_verts_to_choose = num_threads + rand()%(nvertsinC/2);
	}
	#if DEBUG > 1
	printf("Num vertices to choose in L = %d\n", num_verts_to_choose);
	#endif
	while(num_verts_to_choose) {
	    float rnd, cursum = 0.0;
	    rnd = drand48()*sum_of_sel_probability;
	    for (i = 0; i < gr->num_v; i++) {
		if(gr->vl[i].inC && !gr->vl[i].inL) {
		    cursum += 1.0/(2* (gr->vl[i].deg));
		    if (cursum > rnd) {
			gr->vl[i].inL = 1;
			sum_of_sel_probability -= 1.0/(2 * (gr->vl[i].deg));
			#if DEBUG > 1
			printf("Vertex chosen in L = %d\n", i);
			#endif
			break;
		    }
		}
	    }
	    num_verts_to_choose--;
	}

	create_threads(num_threads, mis_transition_from_L_to_S);
	join_threads(num_threads);
#if 0
	#pragma omp parallel for
	for (i = 0; i < gr->num_v; i++) {
	    if(gr->vl[i].inL) {
		if(greatestAmongNeighbors(gr->vl, i)) {
		    /* Put yourself in S. Remove yourself from C */
		    if(gr->vl[i].inS) {
			continue;
		    }
		    gr->vl[i].inS = 1;
		    gr->vl[i].inC = 0;
		    #if DEBUG > 1
		    printf("Putting vertex = %d  weight = %d in S\n", i, gr->vl[i].weight);
		    #endif
		    remVertsPerThread[omp_get_thread_num()]++;
		}
	    }
	}
#endif

	create_threads(num_threads, mis_transition_remove_neighbors_of_S_from_C);
	join_threads(num_threads);
#if 0
	/* Remove neighbors of nodes in S */
#pragma omp parallel for private(j)
	for (i = 0; i < gr->num_v; i++) {
	    if(gr->vl[i].inS) {
		gr->vl[i].inL = 0;
		/* Walk through my neighbor list. Delete each neighbor if I am its greatest neighbor */
		for (j = 0; j < gr->vl[i].neighbor_count; j++) {
		    int idx = gr->vl[i].neighbors[j].id;
		    #if DEBUG > 1
		    printf("vertex = %d neighbor count = %d  Neighbor vertex = %d\n", i, gr->vl[i].neighbor_count, idx);
		    #endif
		    if (gr->vl[idx].inC) {
			if(isX_greatestNeighborOf_Y_inS(gr->vl, i, idx)) {
			    gr->vl[idx].inC = 0;
			    gr->vl[idx].inL = 0;
			    //nvertsinC--;
			    remVertsPerThread[omp_get_thread_num()]++;
			    #if DEBUG > 1
			    printf("vertex = %d  removed vertex = %d from C\n", i, idx);
			    #endif
			}
		    }
		}
	    }
	}
#endif
	for (i = 0; i < num_threads; i++) {
	    sum_removed += thread_data[i];
	}

	nvertsinC -= sum_removed;
    }

    //printf("Total iterations = %d\n", itr);
    //printf("Final set of independent vertices.\n");
    //printinSVerts(gr->vl, gr->num_v);
}


void *do_color_graph(void *tid)
{
    graph *gr = cur_gr;
    int share = (gr->num_v)/num_threads;
    int start = *(int *)(tid) * share;
    int end = start + share - 1;
    int i;
    
    /* last threads gets all remaining */
    if (*(int *)tid == num_threads - 1) {
	end = gr->num_v - 1;
    }

#if DEBUG > 1
    printf("%s: I am thread[%d] working on [%d to %d]\n", __func__, *(int*)tid, start, end);
#endif
    /* no locking needed because each thread
     *  get mutually exclusive share of the vertices
     */
    for (i = start; i <= end; i++) {
	if(gr->vl[i].inS == 1) { /* in MIS */
	    assert(gr->vl[i].color == 0);
	    gr->vl[i].color = get_best_color(i);
	    if (gr->vl[i].color > gr->num_colors) {
		gr->num_colors = gr->vl[i].color;
	    }
	    safe_decrement(&gr->uncolored);
	    gr->vl[i].inS = 0; /* remove from S */
	    gr->vl[i].inC = 0; /* FIXME: needed? remove from C */
	}
    }

    pthread_attr_destroy(&attr);
    pthread_exit(NULL);
}

static void color_graph(graph *gr)
{
    int i;
    while(gr->uncolored > 0) {
	omp_mis_luby(gr);
	/* Apply color for vertices in S*/
	create_threads(num_threads, do_color_graph);
	join_threads(num_threads);
    }

    /* verify that all nodes are colored */
    for (i = 0; i < gr->num_v; i++) {
	assert(gr->vl[i].color > 0);
    }
}

void *do_adjacency_list_next_graph(void *thread_id)
{
    graph *gr = cur_gr;
    graph *next_gr = next_level_gr;
    vertex *vl = gr->vl, *next_vl = next_gr->vl;
    int share = (next_gr->num_v)/num_threads;
    int tid = *(int *)(thread_id);
    int start = tid * share;
    int end = start + share - 1;
    int i, j;

    /* last threads gets all remaining */
    if (tid == num_threads - 1) {
	end = next_gr->num_v - 1;
    }

#if DEBUG > 1
    printf("%s: I am thread[%d] working on [%d to %d]\n", __func__, tid, start, end);
#endif
    /* Build adjacency list for the next graph */
    for (i = start; i <= end; i++) {
	int n = 0;
	next_vl[i].neighbors = calloc(cur_gr->max_degree*2, sizeof(neighbor));
	/* Loop over the predecessors
	 * For each predecessor, loop over its neighbors and add it to the current vertex's neighbor list if not already present
	 */
	for (j = 0; j < 2; j++) {
	    int pred_id = next_vl[i].predecessor_id[j];
	    int k;
	    if (pred_id == -1) continue;
	    
	    for (k = 0; k < vl[pred_id].neighbor_count; k++) {
		int id = vl[pred_id].neighbors[k].id;
		int nl;

		if (vl[pred_id].match == id) continue;

		/* Find if already present in next graph */
		for (nl = 0; nl < next_vl[i].neighbor_count; nl++) {
		    if (next_vl[i].neighbors[nl].id == vl[id].successor_id) {
			next_vl[i].neighbors[nl].weight += vl[pred_id].neighbors[k].weight;
			break;
		    }
		}
		if (nl == next_vl[i].neighbor_count) { /* Not yet present */
		    next_vl[i].neighbors[nl].id = vl[id].successor_id;
		    next_vl[i].neighbors[nl].weight += vl[pred_id].neighbors[k].weight;
		    next_vl[i].neighbor_count++;
		}
	    }
	}
	next_vl[i].deg = next_vl[i].neighbor_count;
    }
    return NULL;
}

/* 
 * A matched edge appears as a vertex in the new graph.
 * Vertex weight is the sum of weights of the two edge vertices.
 * To find neighbors of the new vertex:
 *  1. For each vertex in current graph find its index in the next graph.
 *  2. Those vertices that are mutually matched, map to same index in the next graph.
 *  3. Neighbors of the vertex in the next graph is the union of neighbors of vertices in the prev graph mapping to it.
 *  4. Multiple edges between same vertices in the next graph are mapped to just 1 edge with its weight equal to sum of weights of participating edges.
 */
graph *create_next_level_graph(graph *cur_gr)
{
    graph *next_gr = create_graph();
    int i, j = 0;
    vertex *vl = cur_gr->vl, *next_vl;

    /* number of vertices in next graph = matched verts/2 + unmatched vertices in the current graph */
    next_gr->num_v = (cur_gr->num_v - cur_gr->unmatched)/2 + cur_gr->unmatched;
    next_gr->uncolored = next_gr->num_v;
    next_gr->unmatched = next_gr->num_v; 

    next_gr->vl = (vertex *)calloc(next_gr->num_v, sizeof(vertex));
    next_vl = next_gr->vl;

    /* For each vertex in current graph find its index in the next graph 
     * 1. If vertex is matched and its id is less than that of its match, mark its id in the next graph as j++.
     * 2. If vertex is matched and its id is greater than that of its match, mark its id in the next graph as that of its match in the next graph.
     * 3. If vertex is unmatched, mark its id in the next graph as j++.
     */
    for (i = 0; i < cur_gr->num_v; i++) {
	if(vl[i].match == -1) {/* unmatched */
	    vl[i].successor_id = j;
	    next_vl[j].id = j;
	    next_vl[j].predecessor_id[0] = i;
	    next_vl[j].predecessor_id[1] = -1; /* It has only 1 predecessor */
	    next_vl[j].weight += vl[i].weight;
	    next_vl[j].match = -1;
	    j++; /* should be last */
	}
	else {/* matched */
	    if(i < vl[i].match) { /* id less than that of match */
		vl[i].successor_id = j;
		next_vl[j].id = j;
		next_vl[j].predecessor_id[0] = i;
		next_vl[j].weight += vl[i].weight;
		next_vl[j].match = -1;
		j++; /* should be last */
	    }
	    else {
		int localj = vl[vl[i].match].successor_id;
		vl[i].successor_id = localj;
		next_vl[localj].predecessor_id[1] = i;
		next_vl[localj].weight += vl[i].weight;
	    }
	}
    }
    printf("Current Graph numv = %d, unmatched = %d\n", cur_gr->num_v, cur_gr->unmatched);
    printf("Next Graph j = %d, numv = %d\n", j, next_gr->num_v);
    assert(j == next_gr->num_v);

    next_level_gr = next_gr;
    create_threads(num_threads, do_adjacency_list_next_graph);
    join_threads(num_threads);

    for (i = 0; i < next_gr->num_v; i++) {
	if (next_vl[i].deg > next_gr->max_degree) next_gr->max_degree = next_vl[i].deg;
    }

    cur_gr->next = next_gr;
    return next_gr;
}

long double *x;
void *jacobi(void *thread_id) 
{
    int share = (cur_gr->num_v)/num_threads;
    int tid = *(int *)(thread_id);
    int start = tid * share;
    int end = start + share - 1;
    int i, j;
    
    /* last threads gets all remaining */
    if (tid == num_threads - 1) {
	end = cur_gr->num_v - 1;
    }
    
    for (i = start; i <= end; i++) {
	int temp = 0;
	int Aii = 0;
	for (j = 0; j < cur_gr->vl[i].neighbor_count; j++) {
	    int jj = cur_gr->vl[i].neighbors[j].id;
	    temp -= cur_gr->vl[i].neighbors[j].weight * x[jj];
	    Aii += cur_gr->vl[i].neighbors[j].weight;
	}
	x[i] = -temp /Aii;
    }
    sem_post(&thread_sema);
    return NULL;
}

int main(int argc, char *argv[])
{
    int i, j, it = 0, level = 0;
    struct timeval start_gs, end_gs;
    double time_gs;

    int maxIts = 100000;
    long double *x_old = NULL;
    long double diff;
    long double sum = 0.0;
    long double thresh = 0.0001;
    int zz = 0, nz = 0;    
    if (argc < 2) {
	printf("Usage: ./a.out <graph.txt> <num_threads>\n");
	exit (-1);
    }
    
    if (argc == 3) {
	num_threads = atoi(argv[2]);
	if (num_threads < 1) {
	    printf("specify at least 1 thread \n");
	    exit (-1);
	}
    }

    /* initialize thread ids */
    for (i = 0; i < num_threads; i++) {
      wqueue[i] = wq_init();
      thread_ids[i] = i;
    }
    
    /* initialize graph struct */
    gr_root = create_graph();
    cur_gr = gr_root;
    setup_threads(num_threads);    
    /* populate the graph from given file */
    populate_graph(gr_root, argv[1]);
    
    gettimeofday(&start_gs, NULL);
    
    x_old = (long double *)calloc(cur_gr->num_v, sizeof(long double));
    x = (long double *)malloc(cur_gr->num_v*sizeof(long double));
    
    for (i = 0; i < cur_gr->num_v; i++)
      x[i] = rand() % 10000; //__CHECK__
    
    while (it++ < maxIts)
    {	
#if OMP_PTHREADS
#pragma omp parallel for private(j) num_threads(num_threads)
	for (i = 0; i < cur_gr->num_v; i++) {
	    int temp = 0;
	    int Aii = 0;
	    for (j = 0; j < cur_gr->vl[i].neighbor_count; j++) {
		int jj = cur_gr->vl[i].neighbors[j].id;
		temp -= cur_gr->vl[i].neighbors[j].weight * x[jj]; // __CHECK__
		Aii += cur_gr->vl[i].neighbors[j].weight;
	    }
	    //printf("%d :: %d\n", i, Aii);
	    x[i] = -temp /Aii;
	}
#else
    for (i = 0; i < num_threads; i++) 
      wq_add(wqueue[i], &jacobi);
    for (i = 0; i < num_threads; i++)
      sem_wait(&thread_sema);
#endif	

	sum = 0.0;
	for (i = 0; i < cur_gr->num_v; i++)
	    sum += (x[i] - x_old[i])*(x[i] - x_old[i]);
	
	for (i = 0; i < cur_gr->num_v; i++)
	    x_old[i] = x[i];
	
	sum = sum/cur_gr->num_v;
	diff = sqrt(sum);
//	printf("%Lf, %Lf\n", sum, diff);
	if (diff < thresh)
	    break;
    }
    gettimeofday(&end_gs, NULL);
    time_gs = (end_gs.tv_sec - start_gs.tv_sec);
    printf("GS time = %f\n", time_gs);
    printf("%d iterations to converge with convergence value of %Lf\n", it, diff);
    
//    print_full_graph(cur_gr, 1);
    
#if 0
    while (cur_gr->num_v > 100) {
	#if DEBUG > -1
	printf("======= Level %d ====== \n", level++);
	#endif
	#if DEBUG > 0
	printf("======= Before coloring ====== \n");
	print_full_graph(cur_gr, 0);
	#endif
	
	/* do graph coloring */
	color_graph(cur_gr);

	#if DEBUG > 0
	printf("======= After coloring ====== \n");
	print_full_graph(cur_gr, 0);
	#endif
	
	/* do maximal graph matching */
	maximal_matching(cur_gr);

	#if DEBUG > 0
	printf("======= After matching ====== \n");
	print_full_graph(cur_gr, 1);
	#endif
	
	/* create next level graph */
	cur_gr = create_next_level_graph(cur_gr);
    }

    #if DEBUG > -1
    printf("======= Level %d ====== \n", level);
    printf("======= Final graph ====== \n");
    print_full_graph(cur_gr, 1);
    #endif

#endif

    /* cleanup */
    free_graph(gr_root);
    
    return 0;
}
