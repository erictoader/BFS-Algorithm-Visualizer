# BFS-Algorithm-Visualizer
A command line application that acts as an algorithm visualizer for Breadth-First-Search

! The traversed nodes will point towards their parent nodes using arrows !

## Commands:
- neighb [row] [col] : will display all of the immediate neighbours of a given node
- bfs [row] [col] : will run a BFS traversal starting from the node with the given coordinates
- bfs_step [row] [col] : will run a step-by-step BFS traversal starting from the node with the given coordinates
- bfs_tree [row] [col] : will pretty-print the resulting tree after executing the BFS traversal starting from the given node
- path [row1] [col1] [row2] [col2] : will find and display the shortest path from one node to another using BFS; if the node cannot be reached (i.e. it is blocked by a wall), then the nodes will be coloured in red
- perf : will analyze the complexity of the BFS algorithm, providing you with a detailed analysis as an output .html file
- clear : will reset the grid to its initial blank state
- exit : will terminate the program

## Running bfs_step 6 2
https://user-images.githubusercontent.com/99261319/167446233-8f294125-3840-4aad-959b-cd6ce48f3868.mov

## Running bfs_tree 6 3
<img width="418" alt="Screenshot 2022-05-09 at 18 41 01" src="https://user-images.githubusercontent.com/99261319/167446596-987ffff2-f2fd-4ada-b27d-bdabc975e3c3.png">

## Running path 5 10 3 15
<img width="406" alt="Screenshot 2022-05-09 at 18 42 17" src="https://user-images.githubusercontent.com/99261319/167446761-02117222-55f6-48de-84fd-7e4746b92996.png">
