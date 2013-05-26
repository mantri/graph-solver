#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <assert.h>
#include <pthread.h>
#include <semaphore.h>
#include <float.h>
#include <sys/time.h>
#include "mkl.h"

#define DEBUG -2
#define OMP_PTHREADS 0
#define bool int
#define true 1
#define false 0

// Time Measurements
typedef struct timeval timex;
static timex start, end;
timex checktime() {
	struct timeval now;
	gettimeofday(&now,NULL);
	return now;
}

void starttime() {
	start = checktime();
}

void endtime() {
	end = checktime();
}

long timelapse() {
	long mtime, seconds, useconds;
	endtime();
	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = ((seconds) * 1000 + useconds/1000) + 0.5;
	return mtime;
}

long checktimelapse() {
	long mtime, seconds, useconds;
	timex a = checktime();
	seconds  = a.tv_sec  - a.tv_sec;
	useconds = a.tv_usec - a.tv_usec;
	mtime = ((seconds) * 1000 + useconds/1000) + 0.5;
	return mtime;
}

typedef void*(*funct)(void*);

// Structures 
typedef struct {
	long long int x;
	pthread_mutex_t e_mutex;
	int weight; 
	int valid;
} edge;

typedef struct {
	int id;			/* index of the neighbor */
	int weight;		/* weight of the edge to the  neighbor */
	bool external;
	int root_vid;
} neighbor;

typedef struct {
	int id;				/* integer index */
	int deg;			/* degree of vertex */
	int inC; 
	//int inL;
	int inS;			/* is in the current MIS */
	//int inW;			/* to be considered for next MIS */
	int color; 
	int match; 			/* if matched, its match; else -1 */
	neighbor *neighbors; 		/* adjacency list */
	int neighbor_count; 
	//pthread_mutex_t v_mutex;
	int weight;			
	int predecessor_id[2]; 		/* This vertex's association with atmost 2 vertices in previous graph */
	int successor_id;   		/* This vertex's association with the vertex in next graph */
	int root_vid;
#if DEBUG > 2
	int found;
#endif
} vertex;

typedef struct graph {
	pthread_mutex_t g_mutex;
	int num_v; 			/* number of vertices */
	int num_e; 			/* number of edges */
	vertex *vl; 			/* vertex list */
	edge *edl;  			/* edge list */
	int uncolored; 			/* number of vertices not colored */
	int unmatched; 			/* number of vertices not matched */	
	int num_colors;
	int max_degree;
	struct graph *next;		/* pointer to next graph */
	struct graph *prev;		/* pointer to prev graph */
} graph;

typedef struct queue_elem {
	int v; /* index of the vertex */
	struct queue_elem *next;
} elem;


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


typedef struct queue {
	int size;
	elem *head;
} queue;

// Global Declarations
int root_num_vertices=0;
graph *gr_root;
graph *cur_gr_root;
graph *cur_gr;				/* Points to the current graph being processed. Used by threads for sharing */
graph *next_level_gr;			/* Points to the graph next to current graph being processed. Only used in build_adjacnecy_list_next_graph() */
graph **graph_partitions;
graph *r_next_gr;
queue *q;
queue *qa;
queue *qb;
queue *cur_q;
int *r_max_degree;
int ab;
int sign;
work_queue* wqueue[24];
int num_threads = 8;
int sub_graph_count = 2;
int coarse_size = 100;
//double* Afeast = NULL; // Adjacency Matrix used for spectral bisection
float* A = NULL; // Adjacency Matrix used for spectral bisection
float* I = NULL; // Adjacency Matrix used for spectral bisection
char *bisection = NULL;
pthread_t thread[24];
int thread_ids[24];
int thread_data[24]; // Int Array for thread data. USAGE: memset in your function and access with thread_data[tid].  Create new array for a different data type
pthread_attr_t attr;
pthread_mutex_t q_mutex;
pthread_mutex_t q_mutex_a;
pthread_mutex_t q_mutex_b;
pthread_mutex_t *cur_q_mutex;
int current_color = 0;
float *w;
float *z;
float* fwork;
MKL_INT* iwork;
MKL_INT* ifail;
int should_die = 0;
sem_t mis_sema;
sem_t thread_sema;
struct timeval color_tv1, color_tv2, mis_tv1, mis_tv2, match_tv1, match_tv2, coarsening_tv1, coarsening_tv2;
double color_time[100], mis_time[10000], match_time[100], coarsening_time[100], color_total_time = 0.0, mis_total_time = 0.0, match_total_time = 0.0, coarsening_total_time = 0.0;
int color_itr = 0, mis_itr = 0, match_itr = 0, coarsening_itr = 0;

//QUEUE
static queue* queue_init(pthread_mutex_t *mutex) {
	queue *q = (queue*) malloc(sizeof(queue));
	if (!q) {
		printf("%s: malloc failed \n", __func__);
		exit(-1);
	}
	q->size = 0;
	q->head = NULL;
	pthread_mutex_init(mutex, NULL);
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


static void queue_free(queue *q, pthread_mutex_t *mutex) {
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
	pthread_mutex_destroy(mutex);
}

void wq_free(work_queue *wq) {
	int i;
	assert(wq);
	work *cur = wq->head;
	work *temp;
	while(cur) {
		temp = cur->next;
		free (cur);
		cur = temp;
	}
	free(wq);
	pthread_mutex_destroy(&wq->wq_mutex);
	sem_destroy(&wq->wq_sema);
}


static void queue_add(queue *q, int v, pthread_mutex_t *mutex) {
	elem *e = (elem*) malloc(sizeof(elem));
	if (!e) {
		printf("%s: malloc failed\n", __func__);
		exit(-1);
	}
	e->v = v;
	/* add element into the front of the queue */
	pthread_mutex_lock(mutex);
	e->next = q->head;
	q->head = e;
	q->size++;
	pthread_mutex_unlock(mutex);
}

static int queue_pop(queue *q, pthread_mutex_t *mutex) {
	int ret;
	elem *temp;
	pthread_mutex_lock(mutex);
	ret = q->head->v;
	temp = q->head;
	q->head = q->head->next;
	free(temp);
	q->size--;
	pthread_mutex_unlock(mutex);
	return ret;
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
#if DEBUG > -1
			printf("Thread[%d] exiting...\n", tid);
#endif
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

// PTHREADS
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

// Graph Aux. Functions
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
		if(vl[iy].neighbors[i].external) continue;

		id = vl[iy].neighbors[i].id;
		if (vl[id].inS) {
			if (vl[ix].deg < vl[id].deg) {
#if DEBUG > 3
				printf("vertex = %d not greatest neighbor of vertex = %d in S\n", ix, iy);
#endif
				return 0;
			}
			else if (vl[ix].deg == vl[id].deg) {
				if (vl[ix].weight < vl[id].weight) {
#if DEBUG > 3
					printf("vertex = %d not greatest neighbor of vertex = %d in S\n", ix, iy);
#endif
					return 0;
				}
			}
		}
	}

#if DEBUG > 3
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
		if(vl[id].neighbors[i].external) continue;

		idx = vl[id].neighbors[i].id;
		if (vl[idx].inC) {
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

static void free_graph(graph *gr) {

	assert(gr != gr_root);
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

static void free_graph_levels(graph *gr){
	assert(gr!=gr_root);
	if(gr){
		if(gr->next)
			free_graph_levels(gr->next);
		free_graph(gr);
	}
}

static graph* create_graph() {
	graph *gr = (graph*) malloc(sizeof(graph));
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
	gr->prev = NULL;

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
			printf("%d (%d)", gr->vl[i].neighbors[j].id, gr->vl[i].neighbors[j].external);
		}
		printf("\n");
	}
}

static void populate_graph(graph *gr, char *file) 
{

	FILE *fp;
	int i = 0, x, y;
	int w;
	int j = 0;
	long long e;

	fp = fopen(file, "r");
	fscanf(fp, "%d %d", &gr->num_v, &gr->num_e);
	root_num_vertices = gr->num_v;
	gr->uncolored = gr->num_v;
	gr->unmatched = gr->num_v; 

	// I assumed that there are more vertices than threads. If this assert hits, change the graph or num_threads. The issue will be with sharing across threads.
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
		gr->vl[x].root_vid = x;
		gr->vl[x].deg++;
		gr->vl[x].color = 0; /* colors start from 1 */
		gr->vl[x].match = -1; /* vertex is not yet matched */
		gr->vl[x].inC = 1;
		gr->vl[x].inS = 0;
		//gr->vl[x].inL = 0;
		gr->vl[x].predecessor_id[0] = -1;
		gr->vl[x].predecessor_id[1] = -1;
		gr->vl[x].successor_id = -1;
#if DEBUG > 2
		gr->vl[x].found = -1;
#endif
		if (gr->max_degree < gr->vl[x].deg) { gr->max_degree = gr->vl[x].deg; }

		/* for vertex 'y' */
		gr->vl[y].id = y;
		gr->vl[y].root_vid = y;
		gr->vl[y].deg++;
		gr->vl[y].color = 0;
		gr->vl[y].match = -1;
		gr->vl[y].inC = 1;
		gr->vl[y].inS = 0;
		//gr->vl[y].inL = 0;
		gr->vl[y].predecessor_id[0] = -1;
		gr->vl[y].predecessor_id[1] = -1;
		gr->vl[y].successor_id = -1;
#if DEBUG > 2
		gr->vl[y].found = -1;
#endif
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
		gr->vl[x].neighbors[gr->vl[x].neighbor_count].root_vid = y;
		gr->vl[x].neighbor_count++;

		gr->vl[y].neighbors[gr->vl[y].neighbor_count].id = x;
		gr->vl[y].neighbors[gr->vl[y].neighbor_count].weight = gr->edl[i].weight;	
		gr->vl[y].neighbors[gr->vl[y].neighbor_count].root_vid = x;
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
			if(gr->vl[id].neighbors[j].external) continue;

			idx = gr->vl[id].neighbors[j].id;
			if (gr->vl[idx].color == i) {
				is_seen = 1;
			}
		}
		if (is_seen == 0) 
			return i;
	}
	// printf("%s: ERROR, returning -1 for some reason\n", __func__);
	return i;
}

/* find the neighbor for vertex 'src'
 * that is not yet matched and has the
 * highest weight.
 * New Policy: If my index is greatest of all my neighbors, go ahead else quit.
 * Also find the greatest indexed one of all my neighbors. Then pick that neighbor and see if I am its greatest indexed neighbor.
 * If true pair myself with that neighbor.
 */
int get_best_unmatched_neighbor(graph *gr, int src) {

	int i = 0, ret = -1, max_w = -1, idx = 0;
	for (i = 0; i < gr->vl[src].neighbor_count; i++) {
		if(gr->vl[src].neighbors[i].external) continue;
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
		if(vl[src].neighbors[i].external) continue;
		idx = gr->vl[src].neighbors[i].id;
		if (idx > greatest_neighbor) greatest_neighbor = idx;
	}
	if (src < greatest_neighbor) {
		return -1;	/* let the greatest neighbor handle matching */
	}

	/* Check if I am the greatest neighbor of my greatest neighbor */
	for (i = 0; i < gr->vl[greatest_neighbor].neighbor_count; i++) {
		if(vl[greatest_neighbor].neighbors[i].external) continue;
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
			target = get_best_unmatched_neighbor(gr, src);
			if ( target > -1) {
				gr->vl[src].match = target;
				gr->vl[target].match = src;
				safe_decrement(&gr->unmatched);
				safe_decrement(&gr->unmatched); 
				queue_add(q, src, &q_mutex);
			}
		}
	}

#if DEBUG > 1
	printf("%s: Thread[%d] worked on [%d to %d]\n", __func__, *(int*)tid, start, end);
#endif
	sem_post(&thread_sema);
	return NULL;
}

void maximal_matching(graph *gr)
{
	int colors = gr->num_colors;
	int i;
	elem *cur;
	int src_match, dst_match, dst, j;

	for (i = 1; i <= gr->num_colors; i++) {
		current_color = i;
		q = queue_init(&q_mutex);
#if OMP_PTHREADS == 1
#pragma omp parallel for num_threads(num_threads)
		for (j = 0; j < gr->num_v; j++) {
			int src = j, target;
			/* unmatched vertices of this color */
			if (gr->vl[src].color == i && 
					gr->vl[src].match == -1) {
				target = get_best_unmatched_neighbor(gr, src);
				if (target > -1) {
					gr->vl[src].match = target;
					gr->vl[target].match = src;
					safe_decrement(&gr->unmatched);
					safe_decrement(&gr->unmatched); 
					queue_add(q, src, &q_mutex);
				}
			}
		}

#else
		for (j = 0; j < num_threads; j++) {
#if DEBUG > 1
			printf("Adding do_matching to wqueue[%d]\n", j);
#endif
			wq_add(wqueue[j], &do_matching);
		}

		for (j = 0; j < num_threads; j++)
			sem_wait(&thread_sema);
#endif

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
		queue_free(q, &q_mutex);
	}
}

void * mis_transition_from_C_to_S(void *thread_id)
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
		if(gr->vl[i].inC) {
			if(greatestAmongNeighbors(gr->vl, i)) {
				/* Put yourself in S. Remove yourself from C */
				if(gr->vl[i].inS) {
					continue;
				}
				gr->vl[i].inS = 1;
#if DEBUG > 3
				printf("Putting vertex = %d  weight = %d in S\n", i, gr->vl[i].weight);
#endif
				thread_data[tid]++;
			}
		}
	}

#if DEBUG > 1
	printf("%s: Thread[%d] worked on [%d to %d]\n", __func__, tid, start, end);
#endif
	sem_post(&thread_sema);
	return NULL;
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
			gr->vl[i].inC = 0;
			/* Walk through my neighbor list. Delete each neighbor if I am its greatest neighbor */
			for (j = 0; j < gr->vl[i].neighbor_count; j++) {
				if(gr->vl[i].neighbors[j].external) continue;
				int idx = gr->vl[i].neighbors[j].id;
#if DEBUG > 3
				printf("vertex = %d neighbor count = %d  Neighbor vertex = %d\n", i, gr->vl[i].neighbor_count, idx);
#endif
				if (gr->vl[idx].inC) {
					if(isX_greatestNeighborOf_Y_inS(gr->vl, i, idx)) {
						gr->vl[idx].inC = 0;
						thread_data[tid]++;
#if DEBUG > 3
						printf("vertex = %d  removed vertex = %d from C\n", i, idx);
#endif
					}
				}
			}
		}
	}

#if DEBUG > 1
	printf("%s: Thread[%d] worked on [%d to %d]\n", __func__, tid, start, end);
#endif
	sem_post(&thread_sema);
	return NULL;
}

static void omp_mis_luby(graph *gr)
{
	int i = 0, x, y, nvertsinC = 0, nedges = 0;
	int j = 0, itr = -1;
	long long maxweight, e;
	float sum_of_sel_probability = 0.0;
	int nochange_counter = 0;


	/* Run luby's only for uncolored vertices */
	for (i = 0; i < gr->num_v; i++) {
		if(gr->vl[i].color == 0) { /* not yet colored */
			gr->vl[i].inC = 1; /* include them in the MIS */
			nvertsinC++;
		}
	}

	maxweight = (long long) pow((float) nvertsinC, 4);

	while (nvertsinC) {
		int sum_removed = 0, num_verts_to_choose = 0;
		sum_of_sel_probability = 0.0;
		itr++;

#if DEBUG > 1
		printf("Iteration = %d vertices in C = %d\n", itr, nvertsinC);
#endif

		/*
		   if (nochange_counter > 4) {
		   asm("int $3");  
		   }
		   */

		/* Clear thread_data before usage */
		memset(thread_data, 0, sizeof(thread_data));

		for (i = 0; i < gr->num_v; i++) {
			if(gr->vl[i].inC) {
				gr->vl[i].weight = rand() % maxweight;
				sum_of_sel_probability += 1.0/(2 * (gr->vl[i].deg));
			}
		}

#if OMP_PTHREADS == 1
#pragma omp parallel for num_threads(num_threads)
		for (i = 0; i < gr->num_v; i++) {
			if(gr->vl[i].inC) {
				if(greatestAmongNeighbors(gr->vl, i)) {
					/* Put yourself in S */
					if(gr->vl[i].inS) {
						continue;
					}
					gr->vl[i].inS = 1;
#if DEBUG > 1
					printf("Putting vertex = %d  weight = %d in S\n", i, gr->vl[i].weight);
#endif
					thread_data[omp_get_thread_num()]++;
				}
			}
		}

		/* Remove neighbors of nodes in S */
#pragma omp parallel for private(j) num_threads(num_threads)
		for (i = 0; i < gr->num_v; i++) {
			if(gr->vl[i].inS) {
				gr->vl[i].inC = 0;
				/* Walk through my neighbor list. Delete each neighbor if I am its greatest neighbor */
				for (j = 0; j < gr->vl[i].neighbor_count; j++) {
					if(gr->vl[i].neighbors[j].external) continue;
					int idx = gr->vl[i].neighbors[j].id;
#if DEBUG > 2
					printf("vertex = %d neighbor count = %d  Neighbor vertex = %d\n", i, gr->vl[i].neighbor_count, idx);
#endif
					if (gr->vl[idx].inC) {
						if(isX_greatestNeighborOf_Y_inS(gr->vl, i, idx)) {
							gr->vl[idx].inC = 0;
							thread_data[omp_get_thread_num()]++;
#if DEBUG > 2
							printf("vertex = %d  removed vertex = %d from C\n", i, idx);
#endif
						}
					}
				}
			}
		}

#else
		for (i = 0; i < num_threads; i++) {
#if DEBUG > 1
			printf("Adding mis_transition_from_C_to_S to wqueue[%d]\n", i);
#endif
			wq_add(wqueue[i], &mis_transition_from_C_to_S);
		}

		for (i = 0; i < num_threads; i++) {
			sem_wait(&thread_sema);
		}
		for (i = 0; i < num_threads; i++) {
#if DEBUG > 1
			printf("Adding mis_transition_remove_neighbors_of_S_from_C to wqueue[%d]\n", i);
#endif
			wq_add(wqueue[i], &mis_transition_remove_neighbors_of_S_from_C);
		}

		for (i = 0; i < num_threads; i++) {
			sem_wait(&thread_sema);
		}

#endif
		for (i = 0; i < num_threads; i++) {
			sum_removed += thread_data[i];
		}

		if (sum_removed == 0) {
			nochange_counter++;
		} else {
			nochange_counter = 0;
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
			gr->vl[i].inC = 0; 
		}
	}
#if DEBUG > 1
	printf("%s: Thread[%d] worked on [%d to %d]\n", __func__, *(int*)tid, start, end);
#endif
	sem_post(&thread_sema);
	return NULL;
}

static void color_graph(graph *gr)
{
	int i, sum_colored = 0;

	while(gr->uncolored > 0) {
		gettimeofday(&mis_tv1, NULL);
		omp_mis_luby(gr);
		gettimeofday(&mis_tv2, NULL);
		mis_time[mis_itr++] = (mis_tv2.tv_usec - mis_tv1.tv_usec) + (mis_tv2.tv_sec - mis_tv1.tv_sec)*1000000.0;
		assert(mis_itr < 10000);
		/* Apply color for vertices in S*/
#if OMP_PTHREADS == 1
		memset(thread_data, 0, sizeof(int)*num_threads);
		sum_colored = 0;
#pragma omp parallel for num_threads(num_threads)
		for (i = 0; i < gr->num_v; i++) {
			if(gr->vl[i].inS == 1) { /* in MIS */
				assert(gr->vl[i].color == 0);
				gr->vl[i].color = get_best_color(i);
#if DEBUG > 1
				printf("%s, vertex = %d, color = %d\n", __func__, i, gr->vl[i].color);
#endif
				if (gr->vl[i].color > gr->num_colors) {
					gr->num_colors = gr->vl[i].color;
				}
				thread_data[omp_get_thread_num()]++;
				gr->vl[i].inS = 0; /* remove from S */
				gr->vl[i].inC = 0;
			}
		}
		//#pragma parallel for reduction(+ : sum_colored)
		for (i = 0; i < num_threads; i++) {
			sum_colored += thread_data[i];
		}
		gr->uncolored -= sum_colored;
#else
		for (i = 0; i < num_threads; i++) {
#if DEBUG > 1
			printf("Adding do_color_graph to wqueue[%d]\n", i);
#endif
			wq_add(wqueue[i], &do_color_graph);
		}

		for (i = 0; i < num_threads; i++)
			sem_wait(&thread_sema);

#endif
	}

#if DEBUG >2
	/* verify that all nodes are colored */
	for (i = 0; i < gr->num_v; i++) {
		assert(gr->vl[i].color > 0);
	}
#endif
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
		next_vl[i].neighbors = (neighbor*) calloc(cur_gr->max_degree*2, sizeof(neighbor));
		/* Loop over the predecessors
		 * For each predecessor, loop over its neighbors and add it to the current vertex's neighbor list if not already present
		 */
		for (j = 0; j < 2; j++) {
			int pred_id = next_vl[i].predecessor_id[j];
			int k;
			if (pred_id == -1) continue;

			for (k = 0; k < vl[pred_id].neighbor_count; k++) {
				if(vl[pred_id].neighbors[k].external) continue;
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
					next_vl[i].neighbors[nl].root_vid = vl[id].successor_id;
					next_vl[i].neighbors[nl].weight += vl[pred_id].neighbors[k].weight;
					next_vl[i].neighbor_count++;
				}
			}
		}
		next_vl[i].deg = next_vl[i].neighbor_count;
	}
#if DEBUG > 1
	printf("%s: Thread[%d] worked on [%d to %d]\n", __func__, tid, start, end);
#endif
	sem_post(&thread_sema);
	return NULL;
}

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
		if(vl[i].match == -1) {		/* unmatched */
			vl[i].successor_id = j;
			next_vl[j].id = j;
			next_vl[j].root_vid = vl[i].root_vid;
			next_vl[j].predecessor_id[0] = i;
			next_vl[j].predecessor_id[1] = -1; /* It has only 1 predecessor */
			next_vl[j].weight += vl[i].weight;
			next_vl[j].match = -1;
			j++; /* should be last */
		}
		else {		/* matched */
			if(i < vl[i].match) { /* id less than that of match */
				vl[i].successor_id = j;
				next_vl[j].id = j;
				next_vl[j].root_vid = vl[i].root_vid;
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

#if OMP_PTHREADS == 1
#pragma omp parallel for num_threads(num_threads)
	for (i = 0; i < next_gr->num_v; i++) {
		int n = 0, j;
		next_vl[i].neighbors = calloc(cur_gr->max_degree*2, sizeof(neighbor));
		/* Loop over the predecessors
		 * For each predecessor, loop over its neighbors and add it to the current vertex's neighbor list if not already present
		 */
		for (j = 0; j < 2; j++) {
			int pred_id = next_vl[i].predecessor_id[j];
			int k;
			if (pred_id == -1) continue;

			for (k = 0; k < vl[pred_id].neighbor_count; k++) {
				if(vl[pred_id].neighbors[k].external) continue;
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

#else
	for (i = 0; i < num_threads; i++) {
#if DEBUG > 1
		printf("Adding do_adjacency_list_next_graph to wqueue[%d]\n", i);
#endif
		wq_add(wqueue[i], &do_adjacency_list_next_graph);
	}

	for (i = 0; i < num_threads; i++) {
		sem_wait(&thread_sema);
	}
#endif

	for (i = 0; i < next_gr->num_v; i++) {
		if (next_vl[i].deg > next_gr->max_degree) next_gr->max_degree = next_vl[i].deg;
	}
	cur_gr->next = next_gr;
	next_gr->prev = cur_gr;
	return next_gr;
}

void *populateVertices(void* tid){
	int share = (r_next_gr->num_v)/num_threads;
	int start = *(int *)(tid) * share;
	int end = start + share - 1;
	int i, j, id;
	int *max_degree = &(r_max_degree[*(int*)(tid)]);
	vertex *vl = cur_gr_root->vl; 
	vertex *next_vl = r_next_gr->vl;

	if (*(int *)tid == num_threads - 1) {
		end = r_next_gr->num_v - 1;
	}

	for (i = start; i <= end; i++) {
		id = queue_pop(cur_q, cur_q_mutex); 
		if(r_max_degree[*(int*)(tid)] < vl[id].deg) r_max_degree[*(int*)(tid)] = vl[id].deg;
		vl[id].successor_id = i;
		vl[id].predecessor_id[0] = sign; 
		next_vl[i].id = i;
		next_vl[i].root_vid = vl[id].root_vid;
		next_vl[i].deg = vl[i].deg;
		next_vl[i].inC = 1;
		//next_vl[i].inL = 0;
		next_vl[i].inS = 0;
		next_vl[i].color = 0;
		next_vl[i].match = -1;
		next_vl[i].neighbor_count = vl[id].neighbor_count;
		next_vl[i].neighbors = (neighbor *) calloc(vl[id].neighbor_count, sizeof(neighbor));
		next_vl[i].weight = vl[id].weight;
		next_vl[i].predecessor_id[0] = id;
		next_vl[i].predecessor_id[1] = -1;
		next_vl[i].successor_id = -1;
	}

#if DEBUG > 1
	printf("%s: Thread[%d] worked on [%d to %d]\n", __func__, *(int*)tid, start, end);
#endif
	sem_post(&thread_sema);
	return NULL;
}

void *create_sub_graph_vertices(){
	r_next_gr = create_graph();
	int i, j = 0;
	int id, max_degree=0;
	vertex *vl = cur_gr_root->vl; 
	vertex *next_vl;

	// number of vertices in next graph = matched verts/2 + unmatched vertices in the current graph 
	r_next_gr->num_v = cur_q->size;
	r_next_gr->num_e = -1;
	r_next_gr->vl = (vertex *)calloc(r_next_gr->num_v, sizeof(vertex));
	r_next_gr->uncolored = r_next_gr->num_v;
	r_next_gr->unmatched = r_next_gr->num_v; 
	next_vl = r_next_gr->vl;
	r_next_gr->prev = cur_gr_root;
	r_max_degree = (int*) calloc(num_threads, sizeof(int));

#if OMP_PTHREADS == 1
#pragma omp parallel for private(id) num_threads(num_threads)
	for (i = 0; i < r_next_gr->num_v; i++) {
		id = queue_pop(cur_q, cur_q_mutex); 
		if(r_max_degree[omp_get_thread_num()] < vl[id].deg) r_max_degree[omp_get_thread_num()] = vl[id].deg;
		vl[id].successor_id = i;
		vl[id].predecessor_id[0] = sign; 
		next_vl[i].id = i;
		next_vl[i].root_vid = vl[id].root_vid;
		next_vl[i].deg = vl[i].deg;
		next_vl[i].inC = 1;
		//next_vl[i].inL = 0;
		next_vl[i].inS = 0;
		next_vl[i].color = 0;
		next_vl[i].match = -1;
		next_vl[i].neighbor_count = vl[id].neighbor_count;
		next_vl[i].neighbors = (neighbor *) calloc(vl[id].neighbor_count, sizeof(neighbor));
		next_vl[i].weight = vl[id].weight;
		next_vl[i].predecessor_id[0] = id;
		next_vl[i].predecessor_id[1] = -1;
		next_vl[i].successor_id = -1;
	}
#elif OMP_PTHREADS == 0
	for (j = 0; j < num_threads; j++) wq_add(wqueue[j], &populateVertices);
	for (j = 0; j < num_threads; j++) sem_wait(&thread_sema);

#endif

#if DEBUG > 2
	for (i = 0; i < r_next_gr->num_v; i++) {
		if( next_vl[i].predecessor_id[0] < 0)
			printf("ERROR %s: a predecessor_id was -ve at id = %d\n", __func__, i);
	}
#endif
	for(i = 0; i<num_threads; i++)
		if(max_degree < r_max_degree[i])
			max_degree = r_max_degree[i];
	r_next_gr->max_degree = max_degree;
	free(r_max_degree);
	printf("%s: New polarized graph %d @ = %p\n", __func__, sign, r_next_gr);
	return NULL;
}

void *populateNeighbors(void* tid){
	graph *gr = graph_partitions[ab];
	int share = (gr->num_v)/num_threads;
	int start = *(int *)(tid) * share;
	int end = start + share - 1;
	vertex *pv;
	int i, j;

	if (*(int *)tid == num_threads - 1) {
		end = gr->num_v - 1;
	}

	for(i=start; i <= end; i++){ // For each vertex
		pv = &(gr->prev->vl[gr->vl[i].predecessor_id[0]]);
		for(j = 0; j < gr->vl[i].neighbor_count; j++) { // For each neighbor possible of this vertex
			gr->vl[i].neighbors[j] = pv->neighbors[j]; // Copy the neighbor 
			if( gr->vl[i].neighbors[j].external == false){ // If the neighbor is not external
				gr->vl[i].neighbors[j].external = !( gr->prev->vl[pv->neighbors[j].id].predecessor_id[0] == sign);
				gr->vl[i].neighbors[j].id = gr->prev->vl[pv->neighbors[j].id].successor_id;
			}
#if DEBUG > 1
			printf("%s: %d subgraph vl[%d] neighbor[%d] - external = %d, id = %d, weight = %d\n", __func__, sign, i, j, gr->vl[i].neighbors[j].external, gr->vl[i].neighbors[j].id, gr->vl[i].neighbors[j].weight);
#endif
		}
		gr->vl[i].predecessor_id[0] = -1;
		gr->vl[i].predecessor_id[1] = -1;
	}

#if DEBUG > 1
	printf("%s: Thread[%d] worked on [%d to %d]\n", __func__, *(int*)tid, start, end);
#endif
	sem_post(&thread_sema);
	return NULL;
}


void create_sub_graph_neighbors(){
	int i, j;
	vertex *pv;
	graph *gr = graph_partitions[ab];
	printf("%s: Polarized graph %d  ab = %d @ = %p\n", __func__, sign, ab, gr);
	//for(i=0; i < gr->num_v; i++) if(gr->vl[i].predecessor_id[0] == -1) printf("Found error id = %d\n", i);

#if OMP_PTHREADS == 1
#pragma omp parallel for private(pv, j) num_threads(num_threads)
	for(i=0; i < gr->num_v; i++){ // For each vertex
		pv = &(gr->prev->vl[gr->vl[i].predecessor_id[0]]);
		for(j = 0; j < gr->vl[i].neighbor_count; j++) { // For each neighbor possible of this vertex
			gr->vl[i].neighbors[j] = pv->neighbors[j]; // Copy the neighbor 
			if( gr->vl[i].neighbors[j].external == false){ // If the neighbor is not external
				gr->vl[i].neighbors[j].external = !( gr->prev->vl[pv->neighbors[j].id].predecessor_id[0] == sign);
				gr->vl[i].neighbors[j].id = gr->prev->vl[pv->neighbors[j].id].successor_id;
			}
#if DEBUG > 1
			printf("%s: %d subgraph vl[%d] neighbor[%d] - external = %d, id = %d, weight = %d\n", __func__, sign, i, j, gr->vl[i].neighbors[j].external, gr->vl[i].neighbors[j].id, gr->vl[i].neighbors[j].weight);
#endif
		}
		gr->vl[i].predecessor_id[0] = -1;
		gr->vl[i].predecessor_id[1] = -1;
	}
#elif OMP_PTHREADS == 0
	for (j = 0; j < num_threads; j++) wq_add(wqueue[j], &populateNeighbors);
	for (j = 0; j < num_threads; j++) sem_wait(&thread_sema); 
#endif
}


void add_to_queue(queue *q, pthread_mutex_t *mutex, graph *gr, int vid){
	if(gr == cur_gr_root){
#if DEBUG > 1
		printf("%s: pushing vid = %d\n", __func__, vid);
#endif
		queue_add(q, vid, mutex);
	}
	else {
		if(gr->vl[vid].predecessor_id[0] != -1)
			add_to_queue(q, mutex, gr->prev, gr->vl[vid].predecessor_id[0]);
		if(gr->vl[vid].predecessor_id[1] != -1)
			add_to_queue(q, mutex, gr->prev, gr->vl[vid].predecessor_id[1]);
	}
}


void lapack_ev_solver(){
	int i, j;
	MKL_INT itype = 3;
	char jobz = 'V';
	char range = 'A';
	char uplo = 'U';
	MKL_INT n = cur_gr->num_v;
	//float* a; 
	MKL_INT lda = coarse_size;
	//float* b; 
	MKL_INT ldb = coarse_size;
	float vl = 0.0;
	float vu = 0.0;
	MKL_INT il = 2;
	MKL_INT iu = 2;
	float abstol = 0.01;
	MKL_INT m;
	MKL_INT ldz = coarse_size;
	MKL_INT lwork = 8*cur_gr->num_v;
	MKL_INT info;

	ssygvx(&itype, &jobz, &range, &uplo, &n, A, &lda, I, &ldb, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, fwork, &lwork, iwork, ifail, &info);

#if DEBUG > 2
	printf("Eigen Values\n");
	for(i = 0; i < cur_gr->num_v; i++)
		printf("%f\n", w[i]);
	printf("Eigen Vectors\n");
	for(i = 0; i < cur_gr->num_v; i++){
		for(j = 0; j < cur_gr->num_v; j++){
			printf("%f, ", z[i * coarse_size + j]);
		}
		printf("\n");
	}
#endif

	for (j = 0; j < cur_gr->num_v; j++){
		bisection[j] = ( z[2*coarse_size + j] <= 0) ? -1 : 1;
	}

}


void *constructA(void* tid){
	int share = (cur_gr->num_v)/num_threads;
	int start = *(int *)(tid) * share;
	int end = start + share - 1;
	int i, j, sum;

	if (*(int *)tid == num_threads - 1) {
		end = cur_gr->num_v - 1;
	}

	for (i = start; i <= end; i++) {
		sum = 0;
		for(j = 0; j < cur_gr->vl[i].neighbor_count; j++) {
			if(! cur_gr->vl[i].neighbors[j].external){
				A[coarse_size*i + cur_gr->vl[i].neighbors[j].id] = -cur_gr->vl[i].neighbors[j].weight; 
				sum +=cur_gr->vl[i].neighbors[j].weight;
			}
		}
		A[coarse_size*i + i] = sum;
	}

#if DEBUG > 1
	printf("%s: Thread[%d] worked on [%d to %d]\n", __func__, *(int*)tid, start, end);
#endif
	sem_post(&thread_sema);
	return NULL;
}

void *clearA(void* tid){
	int share = (cur_gr->num_v)/num_threads;
	int start = *(int *)(tid) * share;
	int end = start + share - 1;
	int i, j;

	if (*(int *)tid == num_threads - 1) {
		end = cur_gr->num_v - 1;
	}

	for (i = start; i <= end; i++) 
		for(j = 0; j < cur_gr->num_v; j++)
			A[coarse_size*i + j] = 0;

#if DEBUG > 1
	printf("%s: Thread[%d] worked on [%d to %d]\n", __func__, *(int*)tid, start, end);
#endif
	sem_post(&thread_sema);
	return NULL;
}

void *qbisection(void* tid){
	int share = (cur_gr->num_v)/num_threads;
	int start = *(int *)(tid) * share;
	int end = start + share - 1;
	int i, j;

	if (*(int *)tid == num_threads - 1) {
		end = cur_gr->num_v - 1;
	}

	for(i=start; i<= end; i++){
		if (bisection[i] == 1){
			add_to_queue(qb, &q_mutex_b, cur_gr, i);
		} else if (bisection[i] == -1){
			add_to_queue(qa, &q_mutex_a, cur_gr, i);
		}
	}

#if DEBUG > 1
	printf("%s: Thread[%d] worked on [%d to %d]\n", __func__, *(int*)tid, start, end);
#endif
	sem_post(&thread_sema);
	return NULL;
}

void bisect(int a, int b){
	int i, j;
	void *temp1, *temp2;
	int sum = 0;

#if OMP_PTHREADS == 1
#pragma omp parallel for private(j, sum) num_threads(num_threads)
	for (i = 0; i < cur_gr->num_v; i++) {
		sum = 0;
		for(j = 0; j < cur_gr->vl[i].neighbor_count; j++) {
			if(! cur_gr->vl[i].neighbors[j].external){
				A[coarse_size*i + cur_gr->vl[i].neighbors[j].id] = -cur_gr->vl[i].neighbors[j].weight; 
				sum +=cur_gr->vl[i].neighbors[j].weight;
			}
		}
		A[coarse_size*i + i] = sum;
	}
#elif OMP_PTHREADS == 0
	for (j = 0; j < num_threads; j++) wq_add(wqueue[j], &constructA);
	for (j = 0; j < num_threads; j++) sem_wait(&thread_sema);
#endif

#if DEBUG > 2
	printf("Before Solving for Second EigenVector, Laplacian matrix is:\n");
	for(i=0; i < cur_gr->num_v; i++){
		for(j=0; j < cur_gr->num_v; j++){
			printf("%d\t", (int) A[i*cur_gr->num_v + j]);
		} printf("\n");
	}
#endif

	lapack_ev_solver();

#if OMP_PTHREADS == 1
#pragma omp parallel for private(j) num_threads(num_threads)
	for (i = 0; i < cur_gr->num_v; i++) 
		for(j = 0; j < cur_gr->num_v; j++)
			A[coarse_size*i + j] = 0;


#pragma omp parallel for  num_threads(num_threads)
	for(i=0; i<cur_gr->num_v; i++){
		if (bisection[i] == 1){
			add_to_queue(qb, &q_mutex_b, cur_gr, i);
		} else if (bisection[i] == -1){
			add_to_queue(qa, &q_mutex_a, cur_gr, i);
		}
	}
#elif OMP_PTHREADS == 0
	for (j = 0; j < num_threads; j++) wq_add(wqueue[j], &clearA);
	for (j = 0; j < num_threads; j++) sem_wait(&thread_sema);

	for (j = 0; j < num_threads; j++) wq_add(wqueue[j], &qbisection);
	for (j = 0; j < num_threads; j++) sem_wait(&thread_sema);
#endif


	printf("%s : The new polarized graphs are at indeces -1=%d and 1=%d\n", __func__, a, b);
	temp1 = graph_partitions[a];
	temp2 = graph_partitions[b];
	cur_q = qa; cur_q_mutex = &q_mutex_a; sign = -1;
	create_sub_graph_vertices();
	graph_partitions[a] = r_next_gr;
	cur_q = qb; cur_q_mutex = &q_mutex_b; sign = 1;
	create_sub_graph_vertices();
	graph_partitions[b] = r_next_gr;
	sign = -1; ab = a;
	create_sub_graph_neighbors();
	sign = 1; ab = b;
	create_sub_graph_neighbors();

	if(temp1 && temp1 != gr_root) free_graph(temp1);
	if(temp2 && temp2 != gr_root) free_graph(temp2);
}

void coarsen(){
	int level = 0;
	gettimeofday(&coarsening_tv1, NULL);
	while (cur_gr->num_v > coarse_size) {
		level++;
#if DEBUG > -1
		printf("======= Level %d ====== \n", level);
#endif
#if DEBUG > 0
		printf("======= Before coloring ====== \n");
		print_full_graph(cur_gr, 0);
#endif
		/* do graph coloring */
		gettimeofday(&color_tv1, NULL);
		color_graph(cur_gr); // do graph coloring 
		gettimeofday(&color_tv2, NULL);
		color_time[color_itr++] = (color_tv2.tv_usec - color_tv1.tv_usec) + (color_tv2.tv_sec - color_tv1.tv_sec)*1000000.0;
#if DEBUG > 0
		printf("======= After coloring ====== \n");
		print_full_graph(cur_gr, 0);
#endif
		/* do maximal graph matching */
		gettimeofday(&match_tv1, NULL);
		maximal_matching(cur_gr); // do maximal graph matching 
		gettimeofday(&match_tv2, NULL);
		match_time[match_itr++] = (match_tv2.tv_usec - match_tv1.tv_usec) + (match_tv2.tv_sec - match_tv1.tv_sec)*1000000.0;
#if DEBUG > 0
		printf("======= After matching ====== \n");
		print_full_graph(cur_gr, 1);
#endif
		cur_gr = create_next_level_graph(cur_gr); // create next level graph 
	}
	gettimeofday(&coarsening_tv2, NULL);
	coarsening_time[coarsening_itr++] = (coarsening_tv2.tv_usec - coarsening_tv1.tv_usec) + (coarsening_tv2.tv_sec - coarsening_tv1.tv_sec)*1000000.0;
	printf("======= Final Level Graph %d  Max Degree = %d ====== \n", level, cur_gr->max_degree);
#if DEBUG > -1
	print_full_graph(cur_gr, 1);
#endif
}

void divide_graph(){
	int i, j;
	int num_levels = (int)  log2((double) sub_graph_count) ;
	printf("num_levels = %d\n", num_levels);
	int iter = 0;
	int offset = 0;
	long long time_coarsen = 0; 
	long long time_eigen = 0;
	int a=0, b=0;
	for(i=0; i<num_levels; i++){
		iter = 1 << i;
		offset = sub_graph_count / (1 << (i+1));
		for(j=0; j < iter; j++){
			qa = queue_init(&q_mutex_a);
			qb = queue_init(&q_mutex_b);
			a =  j * (1 << ( num_levels - i));
			cur_gr = graph_partitions[a];
			cur_gr_root = cur_gr;
			b = a + offset;
			printf("%s: Coarsening the graph #%d @%p, num_v = %d\n", __func__, a, cur_gr_root, cur_gr->num_v);
			// Coarsening
			starttime();
			coarsen();
			endtime();
			time_coarsen+= timelapse();
			printf("%s: Coarsening finished for #%d @%p\n", __func__, a, cur_gr_root);

			// Spectral Bisection
			printf("%s: Bisecting the coarsened graph @%p,  a = %d, b = %d\n", __func__, cur_gr, a, b);
			starttime();
			bisect(a, b);
			endtime();
			time_eigen += timelapse();
			printf("%s: Bisecting the coarsened graph @%p,  a = %d @%p, b = %d @%p\n", __func__, cur_gr, a, graph_partitions[a], b, graph_partitions[b]);
			printf("%s: Coarsening and Bisection finished for #%d @%p, - num_v = %d, + num_v = %d\n", __func__, a, cur_gr_root, graph_partitions[a]->num_v, graph_partitions[b]->num_v);

			printf("#########################################################################################################\n");
			// Cleanup
			if(cur_gr_root->next != gr_root)
				free_graph_levels(cur_gr_root->next);
			queue_free(qa, &q_mutex_a);
			queue_free(qb, &q_mutex_b);
		}
	}
	printf("time for coarsening = %ld ms\n", time_coarsen);
	printf("time for bisect = %ld ms\n", time_eigen);
}

void *constructI(void *tid){
	int share = (coarse_size)/num_threads;
	int start = *(int *)(tid) * share;
	int end = start + share - 1;
	int i;

	if (*(int *)tid == num_threads - 1) {
		end = coarse_size - 1;
	}

	for (i = start; i <= end; i++) {
		I[coarse_size*i + i] = 1; 
	}

#if DEBUG > 1
	printf("%s: Thread[%d] worked on [%d to %d]\n", __func__, *(int*)tid, start, end);
#endif
	sem_post(&thread_sema);
	return NULL;
}

long double *x = NULL;
long double *x_old = NULL;
int maxIts = 100000;
int counter = 0;
pthread_mutex_t counter_mutex;
//pthread_mutex_t bar_mutex;
pthread_barrier_t pbar;
unsigned long long bar_ctr = 0;

void *modifiedGS(void *thread_id) 
{
	int tid = *(int *)(thread_id);
	int i, j, t;
	int it = 0; 
	long double diff;
	long double error = 0.0;
	long double thresh = 0.0001;
	bool converged = false;
	for (i = 0; i < graph_partitions[tid]->num_v; i++)
		x[graph_partitions[tid]->vl[i].root_vid] = rand() % 10000;

	while (it++ < maxIts)
	{	
		for (i = 0; i < graph_partitions[tid]->num_v; i++) {
			int temp = 0;
			int Aii = 0;
			for (j = 0; j < graph_partitions[tid]->vl[i].neighbor_count; j++) {
				int jj = graph_partitions[tid]->vl[i].neighbors[j].root_vid;
				temp += graph_partitions[tid]->vl[i].neighbors[j].weight * x[jj];
				Aii += graph_partitions[tid]->vl[i].neighbors[j].weight;
			}
			x[graph_partitions[tid]->vl[i].root_vid] = temp /Aii;
		}
		error = 0.0;
		for (i = 0; i < graph_partitions[tid]->num_v; i++) {
			int vi = graph_partitions[tid]->vl[i].root_vid;
			error += (x[vi] - x_old[vi])*(x[vi] - x_old[vi]);
			x_old[vi] = x[vi];
		}

		error = error/graph_partitions[tid]->num_v;
		diff = sqrt(error);
		//printf("tid: %d, Iteration %d, error = %Lf\n", tid, it, diff);

		if (diff < thresh && converged == false) {
			converged = true;
			pthread_mutex_lock(&counter_mutex);
			counter++;
			pthread_mutex_unlock(&counter_mutex);
		}

		pthread_barrier_wait(&pbar);	

		if (counter == num_threads)
			break;
	}
	printf("%d iterations to converge with convergence value of %Lf\n", it, diff);

	sem_post(&thread_sema);
	return NULL;
}

void solver() {
	int i, j, t;
	int k;
	int tid;
	pthread_mutex_init(&counter_mutex, NULL);    
	//pthread_mutex_init(&bar_mutex, NULL);
	pthread_barrier_init(&pbar, NULL, num_threads);
	x_old = (long double *)malloc(root_num_vertices*sizeof(long double ));
	x = (long double *)calloc(root_num_vertices, sizeof(long double ));

#if OMP_PTHREADS == 1
#pragma omp parallel for private (tid, i, j) num_threads(num_threads)
	for (t = 0; t < num_threads; t++) {
		int it = 0; 
		long double diff;
		long double thresh = 0.0001;
		long double error = 0.0;
		tid = omp_get_thread_num();
		bool converged = false;

		for (i = 0; i < graph_partitions[tid]->num_v; i++)
			x[graph_partitions[tid]->vl[i].root_vid] = rand() % 10000;

		while (it++ < maxIts)
		{	
			for (i = 0; i < graph_partitions[tid]->num_v; i++) {
				int temp = 0;
				int Aii = 0;
				for (j = 0; j < graph_partitions[tid]->vl[i].neighbor_count; j++) {
					int jj = graph_partitions[tid]->vl[i].neighbors[j].root_vid;
					temp += graph_partitions[tid]->vl[i].neighbors[j].weight * x[jj];
					Aii += graph_partitions[tid]->vl[i].neighbors[j].weight;
				}
				x[graph_partitions[tid]->vl[i].root_vid] = temp /Aii;

			}
			error = 0.0;
			for (i = 0; i < graph_partitions[tid]->num_v; i++) {
				int vi = graph_partitions[tid]->vl[i].root_vid;
				error += ((x[vi] - x_old[vi])*(x[vi] - x_old[vi]));
				x_old[vi] = x[vi];
			}
			error = error/graph_partitions[tid]->num_v;
			diff = sqrt(error);

			//printf("tid: %d, Iteration %d, error = %Lf\n", tid, it, diff);

			if (diff < thresh && converged == false) {
				converged = true;
				pthread_mutex_lock(&counter_mutex);
				counter++;
				pthread_mutex_unlock(&counter_mutex);
			}
			pthread_barrier_wait(&pbar);	

			if (counter == num_threads)
				break;

		}
		printf("%d iterations to converge with convergence value of %Lf\n", it, diff); 
	}
#else
	for (i = 0; i < num_threads; i++) 
		wq_add(wqueue[i], &modifiedGS);
	for (i = 0; i < num_threads; i++)
		sem_wait(&thread_sema);
#endif
	//pthread_mutex_destroy(&bar_mutex);
	pthread_mutex_destroy(&counter_mutex);
	pthread_barrier_destroy(&pbar);
}

#if DEBUG > 2
void fine_graph_partition_checker(){
	int i, j, k;
	for(i=0; i<sub_graph_count; i++){
		for(j=0; j<graph_partitions[i]->num_v; j++){
			// Disjoint check
			if(gr_root->vl[graph_partitions[i]->vl[j].root_vid].found != -1){
				printf("%s: ERROR! This vertex was found in some subgraph %d and also in %d\n", __func__, gr_root->vl[graph_partitions[i]->vl[j].root_vid].found, i);
			}
			gr_root->vl[graph_partitions[i]->vl[j].root_vid].found = i;

			// Neighbor count check
			if(gr_root->vl[graph_partitions[i]->vl[j].root_vid].neighbor_count != graph_partitions[i]->vl[j].neighbor_count)
				printf("%s: ERROR! gr_root->vl[%d].neighbor_count = %d, graph_partitions[%d]->vl[%d].neighbor_count = %d\n", __func__, graph_partitions[i]->vl[j].root_vid, gr_root->vl[graph_partitions[i]->vl[j].root_vid].neighbor_count, i, j, graph_partitions[i]->vl[j].neighbor_count);
			for(k=0; k<graph_partitions[i]->vl[j].neighbor_count; k++){
				// Neighbor root vid
				if(graph_partitions[i]->vl[j].neighbors[k].root_vid != gr_root->vl[graph_partitions[i]->vl[j].root_vid].neighbors[k].root_vid){
					printf("%s: ERROR! Neighbor %d of root vertex %d doesn't match graph_partition[%d].vl[%d]\n", __func__, k, graph_partitions[i]->vl[j].neighbors[k].root_vid, i, j);
				}
			}
		}
	}

}
#endif


int main(int argc, char *argv[]) {
	int i, j, k;
	struct timeval start_gs, end_gs;
	double time_gs;
	int nt = 0, omp = 0;

#if OMP_PTHREADS == 1
	omp = 1;
#endif
	if (argc < 4) {
		printf("Usage: ./a.out <graph.txt> <num_threads> <coarse_size> <sub_graphs>\n");
		exit (-1);
	}
	if (argc == 5) {
		num_threads = atoi(argv[2]);
		coarse_size = atoi(argv[3]);
		sub_graph_count = atoi(argv[4]);
		if (num_threads < 1) {
			printf("specify at least 1 thread \n"); exit (-1);
		}
		printf("num_threads = %d, coarse_size = %d omp enabled = %d\n", num_threads, coarse_size, omp);
	}
	// INIT
	gr_root = create_graph();
	graph_partitions = (graph**) calloc(sub_graph_count, sizeof(graph*));
	for (i = 0; i < num_threads; i++) {
		wqueue[i] = wq_init();
		thread_ids[i] = i;
	}
	for (i = 0; i < sub_graph_count; i++) 
		graph_partitions[i] = NULL;

#if OMP_PTHREADS != 1
	setup_threads(num_threads);
#endif
	graph_partitions[0] = gr_root;
	populate_graph(gr_root, argv[1]);
	w = calloc(coarse_size, sizeof(float));
	z = calloc(coarse_size * coarse_size, sizeof(float));
	fwork = calloc(coarse_size * 8, sizeof(float));
	iwork = calloc(5*coarse_size, sizeof(int));
	ifail = calloc(coarse_size, sizeof(int));
	A = (float*) calloc(coarse_size * coarse_size, sizeof(float)); 
	I = (float*) calloc(coarse_size * coarse_size, sizeof(float)); 
	bisection = (char*) calloc(coarse_size, sizeof(char));

#if OMP_PTHREADS == 1
#pragma omp parallel for num_threads(num_threads)
	for (i = 0; i < coarse_size; i++) 
		I[coarse_size*i + i] = 1; 
#elif OMP_PTHREADS == 0
	for (j = 0; j < num_threads; j++) wq_add(wqueue[j], &constructI);
	for (j = 0; j < num_threads; j++) sem_wait(&thread_sema);
#endif

	divide_graph();

#if DEBUG >2
	for (j =0; j<num_threads; j++) 
		nt += graph_partitions[j]->num_v;
	printf("%d, %d\n", nt, root_num_vertices);
	fine_graph_partition_checker();
#endif

	gettimeofday(&start_gs, NULL);
	//solver();
	gettimeofday(&end_gs, NULL);
	time_gs = (end_gs.tv_sec - start_gs.tv_sec);
	printf("modified GS time = %f\n", time_gs);

#if OMP_PTHREADS != 1
	should_die = 1;
	for(i = 0; i < num_threads; i++)
		sem_post(&wqueue[i]->wq_sema);
#endif

	// MLC Timing prints
	for (i = 0; i < mis_itr; i++) mis_total_time += mis_time[i];
	for (i = 0; i < color_itr; i++) color_total_time += color_time[i];
	for (i = 0; i < match_itr; i++) match_total_time += match_time[i];
	for (i = 0; i < coarsening_itr; i++) coarsening_total_time += coarsening_time[i];
	mis_total_time = mis_total_time/1000.0;	
	color_total_time = color_total_time/1000.0;	
	match_total_time = match_total_time/1000.0;	
	coarsening_total_time = coarsening_total_time/1000.0;	
	printf("Coarsening time = %f ms\n", coarsening_total_time);
	printf("MIS total time = %f ms\n", mis_total_time);
	printf("Color total time = %f ms\n", color_total_time);
	printf("Match total time = %f ms\n", match_total_time);

	free(w);
	free(z);
	free(fwork);
	free(iwork);
	free(ifail);
	free(A);
	free(I);
	free(bisection);
	return 0;
}


