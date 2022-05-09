#include <stdlib.h>
#include <string.h>
#include "bfs.h"

Queue* initializeQueue()
{
    Queue* newQueue = (Queue*)malloc(sizeof(Queue));
    newQueue->head = NULL;
    newQueue->tail = NULL;

    return newQueue;
}

QueueNode* initializeNode(Node* toInit)
{
    QueueNode* newNode = (QueueNode*)malloc(sizeof(QueueNode));
    newNode->element = toInit;
    newNode->next = NULL;

    return newNode;
}

void enqueue(Queue* queue, Node* toEnqueue)
{
    QueueNode* enqueuedNode = initializeNode(toEnqueue);

    if (queue->head == NULL)
    {
        queue->head = enqueuedNode;
    }
    else
    {
        queue->tail->next = enqueuedNode;
    }
    queue->tail = enqueuedNode;
}

void dequeue(Queue* queue)
{
    if (queue->head == queue->tail)
    {
        queue->head = NULL;
        queue->tail = NULL;
    }
    else
    {
        QueueNode* dequeuedNode = queue->head;
        queue->head = queue->head->next;

        free(dequeuedNode);
    }
}

int get_neighbors(const Grid *grid, Point p, Point neighb[])
{
    // the point p will have at most 4 neighbors (up, down, left, right)
    // avoid the neighbors that are outside the grid limits or fall into a wall
    // note: the size of the array neighb is guaranteed to be at least 4
    int nrNeighbors = 0;
    Point point;

    if ((p.row - 1 >= 0) && (grid->mat[p.row - 1][p.col] == 0))
    {
        point.row = p.row - 1;
        point.col = p.col;
        neighb[nrNeighbors++] = point;
    }
    if ((p.row + 1 <= grid->rows) && (grid->mat[p.row + 1][p.col] == 0))
    {
        point.row = p.row + 1;
        point.col = p.col;
        neighb[nrNeighbors++] = point;
    }
    if ((p.col - 1 >= 0) && (grid->mat[p.row][p.col - 1] == 0))
    {
        point.row = p.row;
        point.col = p.col - 1;
        neighb[nrNeighbors++] = point;
    }
    if ((p.col + 1 <= grid->cols) && (grid->mat[p.row][p.col + 1] == 0))
    {
        point.row = p.row;
        point.col = p.col + 1;
        neighb[nrNeighbors++] = point;
    }

    return nrNeighbors;
}

void grid_to_graph(const Grid *grid, Graph *graph)
{
    //we need to keep the nodes in a matrix, so we can easily refer to a position in the grid
    Node *nodes[MAX_ROWS][MAX_COLS];
    int i, j, k;
    Point neighb[4];

    //compute how many nodes we have and allocate each node
    graph->nrNodes = 0;
    for(i=0; i<grid->rows; ++i){
        for(j=0; j<grid->cols; ++j){
            if(grid->mat[i][j] == 0){
                nodes[i][j] = (Node*)malloc(sizeof(Node));
                memset(nodes[i][j], 0, sizeof(Node)); //initialize all fields with 0/NULL
                nodes[i][j]->position.row = i;
                nodes[i][j]->position.col = j;
                ++graph->nrNodes;
            }else{
                nodes[i][j] = NULL;
            }
        }
    }
    graph->v = (Node**)malloc(graph->nrNodes * sizeof(Node*));
    k = 0;
    for(i=0; i<grid->rows; ++i){
        for(j=0; j<grid->cols; ++j){
            if(nodes[i][j] != NULL){
                graph->v[k++] = nodes[i][j];
            }
        }
    }

    //compute the adjacency list for each node
    for(i=0; i<graph->nrNodes; ++i){
        graph->v[i]->adjSize = get_neighbors(grid, graph->v[i]->position, neighb);
        if(graph->v[i]->adjSize != 0){
            graph->v[i]->adj = (Node**)malloc(graph->v[i]->adjSize * sizeof(Node*));
            k = 0;
            for(j=0; j<graph->v[i]->adjSize; ++j){
                if( neighb[j].row >= 0 && neighb[j].row < grid->rows &&
                    neighb[j].col >= 0 && neighb[j].col < grid->cols &&
                    grid->mat[neighb[j].row][neighb[j].col] == 0){
                        graph->v[i]->adj[k++] = nodes[neighb[j].row][neighb[j].col];
                }
            }
            if(k < graph->v[i]->adjSize){
                //get_neighbors returned some invalid neighbors
                graph->v[i]->adjSize = k;
                graph->v[i]->adj = (Node**)realloc(graph->v[i]->adj, k * sizeof(Node*));
            }
        }
    }
}

void free_graph(Graph *graph)
{
    if(graph->v != NULL){
        for(int i=0; i<graph->nrNodes; ++i){
            if(graph->v[i] != NULL){
                if(graph->v[i]->adj != NULL){
                    free(graph->v[i]->adj);
                    graph->v[i]->adj = NULL;
                }
                graph->v[i]->adjSize = 0;
                free(graph->v[i]);
                graph->v[i] = NULL;
            }
        }
        free(graph->v);
        graph->v = NULL;
    }
    graph->nrNodes = 0;
}

void bfs(Graph* graph, Node* s, Operation* op)
{
    // at the end of the algorithm, every node reachable from s should have the color BLACK
    // for all the visited nodes, the minimum distance from s (dist) and the parent in the BFS tree should be set
    // for counting the number of operations, the optional op parameter is received
    // since op can be NULL (when we are calling the bfs for display purposes), you should check it before counting:
    // if(op != NULL) op->count();

    Queue* bfsQueue = initializeQueue();
    
    s->color = COLOR_GRAY;
    s->parent = NULL;
    s->dist = 0;

    enqueue(bfsQueue, s);
    if (op != NULL) op->count();

    while (bfsQueue->tail != NULL)
    {
        Node* dequeuedNode = bfsQueue->head->element;
        dequeue(bfsQueue);
        if (op != NULL) op->count();

        for (int i = 0; i < dequeuedNode->adjSize; i++)
        {
            Node* adjacentNode = dequeuedNode->adj[i];
            if (op != NULL) op->count();
            if (adjacentNode->color == COLOR_WHITE)
            {
                enqueue(bfsQueue, adjacentNode);
                adjacentNode->color = COLOR_GRAY;
                adjacentNode->dist = dequeuedNode->dist + 1;
                adjacentNode->parent = dequeuedNode;
                if (op != NULL) op->count(4);
            }
        }

        dequeuedNode->color = COLOR_BLACK;
        if (op != NULL) op->count();
    }

}

void PrettyPrintParentArray(int* parentArray, Point* repr, int numberOfNodes, int nodeIndex, int tabs, bool printHeader)
{
    if (printHeader == true)
    {
        printf("\nPretty-printing the parent array:");
        printHeader = false;
    }
    printf("\n");
    for (int i = 0; i < tabs; i++)
    {
        printf("\t");
    }
    printf("(%d | %d)", repr[nodeIndex].row, repr[nodeIndex].col);
    for (int i = 0; i < numberOfNodes; i++)
    {
        if (parentArray[i] == nodeIndex)
        {
            PrettyPrintParentArray(parentArray, repr, numberOfNodes, i, tabs + 1, printHeader);
        }
    }
}

void print_bfs_tree(Graph *graph)
{
    //first, we will represent the BFS tree as a parent array
    int n = 0; //the number of nodes
    int *p = NULL; //the parent array
    Point *repr = NULL; //the representation for each element in p

    //some of the nodes in graph->v may not have been reached by BFS
    //p and repr will contain only the reachable nodes
    int *transf = (int*)malloc(graph->nrNodes * sizeof(int));
    for(int i=0; i<graph->nrNodes; ++i){
        if(graph->v[i]->color == COLOR_BLACK){
            transf[i] = n;
            ++n;
        }else{
            transf[i] = -1;
        }
    }
    if(n == 0){
        //no BFS tree
        free(transf);
        return;
    }

    int err = 0;
    p = (int*)malloc(n * sizeof(int));
    repr = (Point*)malloc(n * sizeof(Point));
    for(int i=0; i<graph->nrNodes && !err; ++i){
        if(graph->v[i]->color == COLOR_BLACK){
            if(transf[i] < 0 || transf[i] >= n){
                err = 1;
            }else{
                repr[transf[i]] = graph->v[i]->position;
                if(graph->v[i]->parent == NULL){
                    p[transf[i]] = -1;
                }else{
                    err = 1;
                    for(int j=0; j<graph->nrNodes; ++j){
                        if(graph->v[i]->parent == graph->v[j]){
                            if(transf[j] >= 0 && transf[j] < n){
                                p[transf[i]] = transf[j];
                                err = 0;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
    free(transf);
    transf = NULL;

    if(!err){
        // the parrent array is p (p[k] is the parent for node k or -1 if k is the root)
        // when printing the node k, print repr[k] (it contains the row and column for that point)
        // you can adapt the code for transforming and printing multi-way trees from the previous labs
        int root = 0;
        for (int i = 0; i < n; i++)
        {
            if (p[i] == -1)
            {
                root = i;
                break;
            }
        }
        PrettyPrintParentArray(p, repr, n, root, 0, true);
        printf("\n");
    }

    if(p != NULL){
        free(p);
        p = NULL;
    }
    if(repr != NULL){
        free(repr);
        repr = NULL;
    }
}

int shortest_path(Graph *graph, Node *start, Node *end, Node *path[])
{
    bfs(graph, start);

    if (end->parent == NULL)
    {
        return -1;
    }
    else
    {
        int pathLength = 0;
        Node* parentNode = end->parent;
        while (parentNode != NULL)
        {
            path[pathLength++] = parentNode;
            parentNode = parentNode->parent;
        }
        return pathLength;
    }
}

void generateNodes(int vertices, Graph* graph)
{
    for (int i = 0; i < vertices; i++)
    {
        (*graph).v[i]->adjSize = 0;
        (*graph).v[i]->adj = (Node**)malloc(vertices * sizeof(Node));
        for (int j = 0; j < vertices; j++)
        {
            (*graph).v[i]->adj[j] = (Node*)malloc(sizeof(Node));
        }
        (*graph).v[i]->color = COLOR_WHITE;
        (*graph).v[i]->dist = 0;
        (*graph).v[i]->parent = NULL;
    }
}

bool isValidEdge(Node* node1, Node* node2)
{
    if (node1 == node2)
    {
        return false;
    }

    for (int i = 0; i < node1->adjSize; i++)
    {
        if (node1->adj[i] == node2)
        {
            return false;
        }
    }

    return true;
}

void generateEdges(int vertices, int edges, Graph* graph)
{
    srand(time(NULL));

    // making sure the graph is connected
    for (int i = 0; i < vertices - 1; i++)
    {
        (*graph).v[i]->adj[(*graph).v[i]->adjSize] = (*graph).v[i + 1];
        (*graph).v[i + 1]->adj[(*graph).v[i + 1]->adjSize] = (*graph).v[i];
        (*graph).v[i]->adjSize++;
        (*graph).v[i + 1]->adjSize++;
        (*graph).v[i + 1]->parent = (*graph).v[i];
    }

    // randomly creating the remaining edges 
    for (int i = vertices - 1; i < edges; i++)
    {
        int v1 = rand() % vertices;
        int v2 = rand() % vertices;

        if (isValidEdge((*graph).v[v1], (*graph).v[v2]))
        {
            (*graph).v[v1]->adj[(*graph).v[v1]->adjSize] = (*graph).v[v2];
            (*graph).v[v2]->adj[(*graph).v[v2]->adjSize] = (*graph).v[v1];
            (*graph).v[v1]->adjSize++;
            (*graph).v[v2]->adjSize++;
            (*graph).v[v2]->parent = (*graph).v[v1];
        }
        else i--;
    }

}

void performance()
{
    int n, i;
    Profiler p("bfs");

    // vary the number of edges
    for(n=1000; n<=4500; n+=100){
        Operation op = p.createOperation("bfs-edges", n);
        Graph graph;
        graph.nrNodes = 100;
        //initialize the nodes of the graph
        graph.v = (Node**)malloc(graph.nrNodes * sizeof(Node*));
        for(i=0; i<graph.nrNodes; ++i){
            graph.v[i] = (Node*)malloc(sizeof(Node));
            memset(graph.v[i], 0, sizeof(Node));
        }
        generateNodes(100, &graph);
        generateEdges(100, n, &graph);

        bfs(&graph, graph.v[0], &op);
        free_graph(&graph);
    }

    // vary the number of vertices
    for(n=100; n<=200; n+=10){
        Operation op = p.createOperation("bfs-vertices", n);
        Graph graph;
        graph.nrNodes = n;
        //initialize the nodes of the graph
        graph.v = (Node**)malloc(graph.nrNodes * sizeof(Node*));
        for(i=0; i<graph.nrNodes; ++i){
            graph.v[i] = (Node*)malloc(sizeof(Node));
            memset(graph.v[i], 0, sizeof(Node));
        }
        generateNodes(n, &graph);
        generateEdges(n, 4500, &graph);
        // TODO: generate 4500 random edges
        // make sure the generated graph is connected

        bfs(&graph, graph.v[0], &op);
        free_graph(&graph);
    }

    p.showReport();
}
