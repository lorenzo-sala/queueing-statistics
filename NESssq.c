#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define N 10              // max number of machines
#define HETA_ARRIVAL 3000 // Expected value for arrivals
#define HETA_SHORT 40     // Expected value for short station
#define HETA_LONG 960     // Expected value for long station
#define MAX_EVENTS 100000000
#define BETA 0.2 // routing probability
#define MAX_TIME 1000000000000000000

int verbose = 0; 

/*------VARIABLES FOR THE SIMULATION CORE------*/
double sim_clock = 0;
int halt=0;
int event_counter=0;
double heta;
double waiting_short = 0;
double waiting_long = 0;
double beta = BETA;
int nserv_s =0;
int nserv_l =0;
double heta_arrival = 3000;
double heta_short = 40;
double heta_long = 960;

/*------VARIABLES FOR THE REGENERATION METHOD------*/
int cycle_num = 0;
int cycle_in_group = 0;
double total_waiting_short = 0;
double total_waiting_long = 0;
double total_waiting_short_squared = 0;
double S_AA = 0;
double S_Anu = 0;
double S_nu = 0;
double S_nunu = 0;
double error = 0;
double r_hat = 0;

typedef enum {
    AS = 1,
    AL = 2,
    DS = 3,
    DL = 4,
    END = 5
} event_label;
/* Definition of the type used to specify the pointer to a node of a list or queue */
typedef struct node* nodePtr;
/*dobbiamo inizializzarlo prima prechè
poi quando definiamo il nodo dobbiamo usare 
nodePtr*/

typedef struct DLL{
    nodePtr Head;
    nodePtr Tail;
} dll;


dll FEL = {NULL, NULL}; /* Pointer to the header of the Future Event List implemented with a doubly linked list */
dll IQ1 = {NULL, NULL}; /* Pointer to the header of the Future Event List implemented with a doubly linked list */
dll IQ2 = {NULL, NULL};


/* Definition of the Event Notice - typical fields to contain event’s attributes */
typedef struct {
    event_label type;
    int machine_id;
    double create_time;
    double occur_time; 
    double arrival_time;  
    double service_time; 
} event_notice;

/* Definition of the Node for managing Event Notices in Lists and Queues */
struct node {
    event_notice event;    // Contains the event data
    nodePtr left;          // Pointer to the previous node in the doubly linked list
    nodePtr right;         // Pointer to the next node in the doubly linked list
};

/* function declaration */
nodePtr get_new_node();
double exponential_random(double heta);
void initialize();
void engine(void);
void a_short(struct node* node_event);
void a_long(struct node* node_event);
void d_short(struct node* node_event);
void d_long(struct node* node_event);
void end(struct node* node_event);
void release_nodes(nodePtr *head);
void report(double sim_duration);
void schedule(struct node* node_event);
struct node* event_pop(void);
void enqueue(dll* IQ, struct node* new_inqueue);
struct node* dequeue(dll* IQ);
int RegPoint(nodePtr node_event, int *cycle_in_group, int *cycle_num);
void CollectRegStatistics(
    double *waiting_long, 
    double *total_waiting_long, 
    double *S_AA,
    double *S_Anu,
    double *S_nu,
    double *S_nunu,
    int *cycle_in_group
    );
void ResetMeasures(double *waiting_long, int *cycle_in_group);
double ComputeConfidenceIntervals(
    double *r_hat, 
    double total_waiting_long,
    int cycle_num, 
    double S_AA, 
    double S_Anu, 
    double S_nunu, 
    double S_nu
);
int DecideToStop(int cycle_num, double error, double r_hat);



nodePtr get_new_node() {
    // Allocate memory for the new node
    nodePtr new_node = (nodePtr)malloc(sizeof(struct node));

    // Check if memory allocation was successful
    if (new_node == NULL) {
        printf("Error: Memory allocation failed!\n");
        exit(EXIT_FAILURE); // Exit the program if malloc fails
    }

    // Initialize the new node with the given event
    new_node->event.type = 0;             // Default type
    new_node->event.machine_id = 0;    // Empty int
    new_node->event.create_time = 0.0;   // Default time
    new_node->event.occur_time = 0.0;    // Default time
    new_node->event.arrival_time = 0.0;  // Default time
    new_node->event.service_time = 0.0;  // Default time
    new_node->left = NULL;  // Set the left pointer to NULL
    new_node->right = NULL; // Set the right pointer to NULL

    return new_node; // Return the pointer to the new node
};

double exponential_random(double heta) {
  double unif = (double)rand() / RAND_MAX;
  if ((1 - unif) < 1e-20) {
    unif = 1e-20; // this will give us a limit in the case in which the argument
                  // of the logarithm is too close to 0 (and thus exploding)
  }
  double exp = (-heta * logf(1 - unif));
  // double pick = 1- ((double)rand() / RAND_MAX);
  return exp;
}

/*
-----------------------------------------------------------------------------
-------------------------SIMULATION CORE-------------------------------------
-----------------------------------------------------------------------------
*/

/*function to print fel for debugging*/
void print_fel() {
    struct node* temp = FEL.Head;
    printf("FEL: ");
    while (temp != NULL) {
        printf("[Event Time: %f, Machine ID: %d, Event Type: %d] -> ", temp->event.occur_time, temp->event.machine_id, temp->event.type);
        temp = temp->right;
    }
    printf("NULL\n\n");
}

void initialize() {
    /* Definition and Initialization of Useful Variables */
    sim_clock = 0;
    halt = 0;

    /*we pick the 10 arrival times of the system. we create the nodes and schedule them!*/
    for (int i = 1; i <= 10; i++) {
        nodePtr init_job = get_new_node();
        double arrival_t = exponential_random(heta_arrival);

        init_job->event.arrival_time = arrival_t;
        init_job->event.occur_time = arrival_t;
        init_job->event.create_time = sim_clock;
        init_job->event.type = AS;
        init_job->event.machine_id = i;
        printf("arrival time picked as %f for machine %d\n",init_job->event.occur_time, init_job->event.machine_id);
        schedule(init_job); 
    }

    verbose ? print_fel() : 0;
   
    /* Initialize Event notice of End Simulation and Schedule last event */
    nodePtr end_event_notice = get_new_node(); 
    end_event_notice->event.type = END;  // Set the type to END
    end_event_notice->event.occur_time = MAX_TIME;  // Set the end time (assuming you have this predefined)
    schedule(end_event_notice);  // Schedule the end simulation event in the Future Event List
}

void engine(void){
    char event_type;
    double old_sim_clock;
    double interval;
    nodePtr new_event; // pointer to struct node

    /* get the 1st event from future event list */
    new_event = event_pop(); 

    /* update clock */
    old_sim_clock = sim_clock;
    sim_clock = new_event->event.occur_time;
    interval = sim_clock - old_sim_clock;

    verbose ? printf("Now processing an event of type %d with occur time %f for machine %d\n", new_event->event.type, new_event->event.occur_time, new_event->event.machine_id) : 0;

    if (RegPoint(new_event, &cycle_in_group, &cycle_num)) {
        verbose ? printf("Found a regeneration point. Ended cycle number %d\n", cycle_num) : 0;
        CollectRegStatistics(
            &waiting_long, 
            &total_waiting_long, 
            &S_AA, 
            &S_Anu, 
            &S_nu, 
            &S_nunu, 
            &cycle_in_group
            );
        ResetMeasures(
            &waiting_long, 
            &cycle_in_group
        );
        error = ComputeConfidenceIntervals(
            &r_hat, 
            total_waiting_long, 
            cycle_num, 
            S_AA, 
            S_Anu, 
            S_nunu, 
            S_nu
            );
    }

    /* identify and process current event */
    event_type = new_event->event.type;
    switch(event_type) {
        case AS : 
        verbose ? printf("Now processing an arrival to the short station (event %d)\n", new_event->event.type) : 0;
        a_short(new_event);
        break;
        case AL : 
        verbose ? printf("Now processing an arrival to the long station (event %d)\n", new_event->event.type) : 0;
        a_long(new_event);
        break;
        case DS : 
        verbose ? printf("Now processing a departure from the short station (event %d)\n", new_event->event.type) : 0;
        d_short(new_event);
        break;
        case DL : 
        //verbose ? printf("Now processing a departure from the station (event %d)\n", new_event->event.type) : 0;
        d_long(new_event);
        break;
        case END : 
        verbose ? printf("Now processing the end (event %d)\n", new_event->event.type) : 0;
        end(new_event);
        break;
    }
    verbose ? print_fel() : 0;
    verbose ? printf("Number of clients in Queue 1: %d\nNumber of clients in Queue 2: %d\n\n", nserv_s, nserv_l) : 0;
    event_counter++;
}

void a_short(struct node* node_event){
    struct node* next_job;
    nserv_s++;
    node_event->event.create_time = sim_clock;

    /*pick a service time*/
    node_event->event.service_time = exponential_random(heta_short);
    printf("Extracted a short departure time of %f\n", node_event->event.service_time);

    /*idle... or busy?*/
    if (nserv_s == 1) { //idle!
        node_event->event.type = DS;
        node_event->event.occur_time = sim_clock + node_event->event.service_time;
        schedule(node_event);
    } else { //busy!
        node_event->event.type = DS;
        node_event->event.occur_time = 0; //initialize it
        enqueue(&IQ1, node_event);
    }
}

void a_long(struct node* node_event) {
    struct node* next_job;
    nserv_l++;
    node_event->event.create_time = sim_clock;

    /*pick a service time*/
    node_event->event.service_time = exponential_random(heta_long);
    printf("Extracted a long departure time of %f\n", node_event->event.service_time);

    /*idle... or busy?*/
    if (nserv_l == 1) { //idle!
        node_event->event.type = DL;
        node_event->event.occur_time = sim_clock + node_event->event.service_time;
        schedule(node_event);
    } else { //busy!
        node_event->event.type = DL;
        node_event->event.occur_time = 0; //initialize it
        enqueue(&IQ2, node_event);
    }
}

void d_short(struct node* node_event){
    double waiting_time;
    nodePtr next_job;

    /*compute our statistics and waiting times!*/
    nserv_s--;
    waiting_time = sim_clock - node_event->event.arrival_time;
    waiting_short = waiting_short + waiting_time; 

    /*decide whether the customer will go to long or short repair station*/
    double decide_route = (rand() / (double)RAND_MAX);
    if (decide_route <= beta) {
        /* go to the long station.
        this means that our d_short becomes a a_long*/
        node_event->event.type = AL;
        node_event->event.arrival_time = node_event->event.occur_time;
        verbose ? printf("Bad luck! Machine %d gets a long repair!\n", node_event->event.machine_id) : 0;
        schedule(node_event);
        if (nserv_s>0) {
            verbose ? printf("queue was empty. dequeuing next job from the queue 1\n") : 0;
            /*extract from queue of server 1 and schedule departure*/
            next_job = dequeue(&IQ1);
            next_job->event.type = DS;
            /*pick service time*/
            next_job->event.service_time = exponential_random(heta_short);
            next_job->event.occur_time = sim_clock + next_job->event.service_time;
            schedule(next_job);
        }

    } else {
        verbose ? printf("Good luck! Machine %d goes back in the pool!\n", node_event->event.machine_id) : 0;
        /*go back to the pool. schedule its next arrival at short station*/
        node_event->event.type = AS;
        node_event->event.arrival_time = sim_clock + exponential_random(heta_arrival);
        node_event->event.occur_time = node_event->event.arrival_time;
        schedule(node_event);
    }
}

void d_long(struct node* node_event){
    double waiting_time;
    struct node* next_job;

    /*compute our statistics and waiting times!*/
    nserv_l--;
    waiting_time = sim_clock - node_event->event.arrival_time;
    waiting_long = waiting_long + waiting_time; 

    node_event = dequeue(&IQ2);
    /*go back to the pool. schedule its next arrival at short station*/
    node_event->event.type = AS;
    double arrival_time = exponential_random(heta_arrival);
    printf("The machine %d goes back to the station. It will break in %f time units\n", node_event->event.machine_id, arrival_time);
    node_event->event.arrival_time = sim_clock + arrival_time;
    node_event->event.occur_time = node_event->event.arrival_time;
    schedule(node_event);
}

void end(struct node* node_event){
    halt = 1;
}

void release_nodes(nodePtr *head) {
    /*traverse the lists and free the nodes*/
    nodePtr temp;
    while (*head != NULL) {
        temp = *head;
        *head = (*head)->right;
        free(temp);
    }
}

void report(double sim_duration) {
    printf("\n==========================================\n");
    printf("Simulation complete!");
    printf("\n==========================================\n");
    printf("\nExecution time: %f seconds\n", sim_duration);
    printf("Number of events processed: %d\n", event_counter);
    printf("Number of regeneration cycles: %d", cycle_num);
    printf("Observation period: %f time units\n", sim_clock);
    printf("Average waiting time at the long repair station: %f\n", r_hat);
    printf("Confidence interval at 0.95 level: (%f, %f)\n", r_hat-error, r_hat+error);
    printf("The error is %f %% of the value of the average", 2*error/r_hat);

    release_nodes(&FEL.Head);
    release_nodes(&IQ1.Head);
    release_nodes(&IQ2.Head);
}

/*
|-------------------------------------------------------------------------------|
|--------------------------------DATA STRUCTURES--------------------------------|
|-------------------------------------------------------------------------------|
*/


void schedule(struct node* node_event){
/* empty fel, head is null */
if (FEL.Head == NULL) {
    FEL.Head = node_event;
    FEL.Tail = node_event;
    node_event->left = NULL;
    node_event->right = NULL;
    verbose ? printf("Scheduled event of type %d at occur time %f for machine %d. EMPTY Q, INSERT HEAD\n", node_event->event.type, node_event->event.occur_time, node_event->event.machine_id) : 0;
    return;
}

/*insert at the front*/
if (node_event->event.occur_time <= FEL.Head->event.occur_time) {
    node_event->right = FEL.Head;
    node_event->left = NULL;
    FEL.Head->left = node_event;
    FEL.Head = node_event;
    verbose ? printf("Scheduled event of type %d at occur time %f for machine %d. INSERTION AT THE FRONT\n", node_event->event.type, node_event->event.occur_time, node_event->event.machine_id) : 0;
    return;
}

/*insert in the middle or end: needed traversal*/

struct node* current = FEL.Head;

/*traverse*/
while (current != NULL && node_event->event.occur_time > current->event.occur_time) {
    current = current->right;
}

/*insert*/
if (current != NULL) {
    node_event->right = current;
    node_event->left = current->left;
    if (current->left != NULL) {
        current->left->right = node_event;
    }
    current->left = node_event;
} else {
    /*insert at the end*/
    node_event->left = FEL.Tail;
    node_event->right = NULL;
    FEL.Tail->right = node_event;
    FEL.Tail = node_event;
}
verbose ? printf("Scheduled event of type %d at occur time %f for machine %d\n", node_event->event.type, node_event->event.occur_time, node_event->event.machine_id) : 0;
}

struct node* event_pop(void){
/*check if the FEL is empty*/
if (FEL.Head == NULL) {
    printf("the future event list is empty. I popped all there was to pop!");
    return NULL;
}

/*save temporarily first node*/
struct node* first_node = FEL.Head;

/*update the head of the FEL to point to the next node*/
FEL.Head = FEL.Head->right;

/*if the list becomes empty update the tail of the FEL*/
if (FEL.Head == NULL) {
    FEL.Tail = NULL; // no nodes left
} else {
    FEL.Head->left = NULL; //detach popped node
}

/*detach pointers*/
first_node->left = NULL;
first_node->right = NULL;

verbose ? printf("Popped node from the FEL. This event was of the type %d and had an occur time of %f\n", first_node->event.type, first_node->event.occur_time) : 0;
return first_node;
}

void enqueue(dll* IQ, struct node* new_inqueue){
    /*check if the input node is null*/
    if (new_inqueue == NULL) {
        printf("why did you give me a null node?");
        return;
    }

    /*init pointers of the new_inqueue*/
    new_inqueue->left = NULL;
    new_inqueue->right = NULL;

    if (IQ->Head == NULL && IQ->Tail == NULL) {
        /*if the queue is empty insert first node*/
        IQ->Head = new_inqueue;
        IQ->Tail = new_inqueue;
        verbose ? printf("The queue was empty. Inserted this node about maching %d at the front of the queue.\n", new_inqueue->event.machine_id) : 0;
    } else {
        /*link new node at the end of the Q*/
        new_inqueue->left = IQ->Tail; // its previous node is the tail
        IQ->Tail->right = new_inqueue; // it is the successor of the tail
        IQ->Tail = new_inqueue; // it is the new tail!
        verbose ? printf("The queue was not empty. Inserted this node about maching %d at the end of the queue.\n", new_inqueue->event.machine_id) : 0;
    }
}

struct node* dequeue(dll* IQ) {
    /*check if valid q*/
   
   /*
    if (IQ == NULL && IQ->Head == NULL) {
        printf("I can't dequeue from an empty/nonexistent queue!");
        return NULL;
    }*/

    /*get first node*/
    struct node* dequeued_node = IQ->Head;
   

    /*update head to point to the second node*/
    IQ->Head = dequeued_node->right;
     printf("reached dequeue function\n");
    

    if (IQ->Head != NULL) {
        /*update the head if there are events left in the q*/
        IQ->Head->left = NULL;
    } else {
        /*if the queue is empty also the tail is null*/
        IQ->Tail = NULL;
    }

    /*disconnect the dequeued node from the queue*/
    dequeued_node->left = NULL;
    dequeued_node->right = NULL;

    printf("dequeued node about machine %d.\n", dequeued_node->event.machine_id);
    return(dequeued_node);
}

/*
|--------------------------------------------------------------------------|
|--------------------------REGENERATION UTILITIES--------------------------|
|--------------------------------------------------------------------------|
*/

int RegPoint(nodePtr node_event, int *cycle_in_group, int *cycle_num) {
    /*
    in this scenario every departure from any station may be a suitable regeneration point. 
    since we are interested in the departure from the long station we pick as regeneration points the departures from the long station.
    in order to preserve the conditions to apply the central limit theorem we have to have a reasonable sample size and so we must group different regeneration cycles together.
    we choose 50 as our sample size, since the minimal number of samples commonly used as guideline is 30.

    this function returns 1 (true) for regeneration point, 0 (false) otherwise
    */
   if (*cycle_in_group < 50)
   {
    if (node_event->event.type == DL) {
    (*cycle_in_group)++;
    return 0;
    }
   } else {
    (*cycle_num)++;
    return 1;
   }
   return 0;
}

void CollectRegStatistics(
    double *waiting_long, 
    double *total_waiting_long, 
    double *S_AA,
    double *S_Anu,
    double *S_nu,
    double *S_nunu,
    int *cycle_in_group
    ) {
    *total_waiting_long += *waiting_long;
    *S_AA += pow(*waiting_long, 2);
    *S_Anu += *waiting_long * *cycle_in_group;
    *S_nu += *cycle_in_group;
    *S_nunu += pow(*cycle_in_group, 2);
}

void ResetMeasures(double *waiting_long, int *cycle_in_group){
    *waiting_long = 0;
    *cycle_in_group = 0;
}

double ComputeConfidenceIntervals(
    double *r_hat, 
    double total_waiting_long,
    int cycle_num, 
    double S_AA, 
    double S_Anu, 
    double S_nunu, 
    double S_nu
    ) {
    *r_hat = total_waiting_long / S_nu;
    double delta = sqrt(cycle_num/(cycle_num-1))*(sqrt(S_AA-2* *r_hat* S_Anu+pow(*r_hat, 2)* S_nunu))/(S_nu);
    error = 1.96 * delta;
    return(error);
}

int DecideToStop(int cycle_num, double error, double r_hat){
    double error_percentage = 2*error/r_hat;
    return(cycle_num > 40 && error_percentage < 0.10);  
}

int main(int argc, char *argv[]){
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] == '-' && argv[i][1] == 'v') {
        verbose = 1; // Enable verbose logging
        break;       // No need to check further once we find the flag
        }
    }
    srand(time(NULL)); 
    initialize(); // do the initialization
    printf("\n\nfinished initialization\n\n");
    clock_t start_time = clock(); // start the stopwatch
    
    /*simulate*/
    while (halt == 0 //&& DecideToStop(cycle_num, error, r_hat) == 0
    ) {
        engine();
    }
    
    clock_t end_time = clock();
    // Calculate the time taken and print it
    double sim_duration = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    report(sim_duration);

    return 0;
}