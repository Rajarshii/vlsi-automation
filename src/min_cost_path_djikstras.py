import heapq
import sys
#import pprint
import numpy as np

# Define directions 
dx = [-1, 1, 0, 0]
dy = [0, 0, -1, 1]

# Check if within the grid bounds
def is_valid_idx(x,y,rows,cols):
    return 0 <= x < rows and 0 <= y < cols

# Check if neighbor
def is_neighbor_bool(x,y):
    #neighbours = [(nx,ny) for (nx,ny) in (x+dx,y+dy)]
    neighbors = []
    for i in range(4):
        neighbors.append((x[0]+dx[i],x[1]+dy[i]))
    return y in neighbors

# Return the neighbors of a cell - we use it only for the backtracking
# the argument group3 is named as such to indicate that
def get_neighbors (x,y,group3):
    neighbors = []
    for i in range(4):
        neighbors.append((x+dx[i],y+dy[i]))  
    #print(neighbors) 
    matching_neighbors = [neighbor for neighbor in neighbors if neighbor in group3]
    #print(matching_neighbors)
    return matching_neighbors 

# Find the minim cost neighbor
# Exclude the neighbor if it is already in path
def get_min_neighbor (neighbors,cost,path):
    min_cost = np.Infinity
    for neighbor in neighbors:
        x_temp,y_temp = neighbor
        if cost[x_temp][y_temp] <= min_cost :
            min_cost = cost[x_temp][y_temp]
            if (x_temp, y_temp) in path:
                continue
            x,y = x_temp,y_temp
    return x,y

# Check if cell inside free island
def is_inside_free_island(coordinate, free_island):
    x, y = coordinate
    x1, y1, x2, y2 = free_island
    return x1 <= x <= x2 and y1 <= y <= y2

# Main algorithm implementation - run Dijkstra's algorithm for grids, and find the minimum cost paths
def dijkstra(grid, source, target, rows, cols, island):
    #rows, cols = len(grid), len(grid[0])
    min_cost = [[float('inf')] * cols for _ in range(rows)]
    #print(source)
    source_row, source_col = source
    min_cost[source_row][source_col] = 0
    target_row, target_col = target

    pq = [(0, source)]
    g3 = []
    while pq:
        cost, (x, y) = heapq.heappop(pq)
        g3.append((x,y))
        if (x, y) == (target_row,target_col):
            break

        if cost > min_cost[x][y]:
            continue

        for i in range(4):
            nx, ny = x + dx[i], y + dy[i]

            if 0 <= nx < rows and 0 <= ny < cols and grid[nx][ny] != np.Infinity:
                new_cost = min_cost[x][y] + 1
                #if grid[nx][ny] == 0:
                if is_inside_free_island((nx,ny),island) and is_inside_free_island((x,y),island):
                    new_cost = min_cost[x][y]
                #elif not is_inside_free_island((x,y),island) and is_inside_free_island((x,y),island):
                #    new_cost = min_cost[x][y] + 1

                if new_cost < min_cost[nx][ny]:
                    min_cost[nx][ny] = new_cost
                    heapq.heappush(pq, (new_cost, (nx, ny)))

    return min_cost, g3


def find_path(grid, source, target, island):
    rows, cols = len(grid), len(grid[0])
    g3 = []
    min_cost, g3 = dijkstra(grid, source, target, rows, cols, island)
    if min_cost[target[0]][target[1]] == float('inf'):
        return np.Infinity, "No path exists."

    # Initialize the list for path, and add target to path
    path = []
    x, y = target
    g3.reverse()
    (x,y) = g3[0]
    path.append((x,y))

    while (x, y) != source:
        neighbors = get_neighbors(x,y,g3)
        x,y = get_min_neighbor(neighbors,min_cost,path)
        #print((x,y))
        path.append((x,y))
    path.reverse()

    # We were working with decremented values, update it here before returning
    for cell in range(len(path)):
        path[cell] = np.add(path[cell],(1,1))

    min_cost_to_target = min_cost[target[0]][target[1]]
    return min_cost_to_target,path



def main(filename):
    with open(filename, "r") as file:
        rows, cols = map(int, file.readline().split())
        source = tuple(np.subtract(tuple(map(int, file.readline().split())),(1,1)))
        target = tuple(np.subtract(tuple(map(int, file.readline().split())),(1,1)))
        blockage = tuple(map(int, file.readline().split()))
        free_island = tuple(np.subtract(tuple(map(int, file.readline().split())),(1,1,1,1)))
    
    grid = [[1] * cols for _ in range(rows)]
    for i in range(blockage[0] - 1, blockage[2]):
        for j in range(blockage[1] - 1, blockage[3]):
            grid[i][j] = float('inf')
    
    for i in range(free_island[0], free_island[2]+1):
        for j in range(free_island[1], free_island[3]+1):
            grid[i][j] = 0
            
    #pprint.pprint(grid)
    total_cost, path = find_path(grid, source, target, free_island)
    if isinstance(path, str):
        print(path)
    else:
        #total_cost = len(path) - 1
        print(f"Minimum cost: {total_cost}")
        print("Path:")
        for x, y in path:
            print(f"({x}, {y})",end="")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 das2541.py input_filename")
        sys.exit(1)
    input_filename = sys.argv[1]
    main(input_filename)
