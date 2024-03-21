import sys
#import numpy as np
import time


def read_file(filename) :
    with open(filename, 'r') as file:
        # note: delete whitespace characters from both ends before reading
        num_nodes, num_edges = map(int, ((file.readline()).strip()).split())
        edges = []  # define the edge as a list
        nodes = {}  # Used for extra check - might not be needed
        graph_nodes = [] # Used to store the actual nodes passed to the graph class

        #debug#print(f"Number of nodes = {num_nodes} | Number of edges = {num_edges}")
        # Iterate over each line in file and split the file lines to get u and v
        #for node in range(num_nodes) :

        for line in file:
            u, v = map(int, (line.strip()).split())
            if u not in nodes: # This can be removed, retained for safe-keeping
                nodes[u] = node(u)
            if v not in nodes:
                nodes[v] = node(v)

            edges.append(edge(u, v))
            #print(f"Added edge : ({u},{v}))")
        for node_ in range(1,num_nodes+1) :
            graph_nodes.append(node(node_)) 

        #debug#print(f"size of edge = {len(edges)} , size of nodes : {len(nodes)}, size of graph_nodes : {num_nodes}")

        if(len(nodes) != len(graph_nodes)) :
            print(f"ERROR: The number of nodes in line 1 : {num_nodes} and nodes in file : {len(nodes)} are different")
            print("ERROR: The nodes are expected to be continuous. Exiting...")
            sys.exit(1)  # exit teh system , erroneous condition
        
        #return graph(list(nodes.values()), edges)
        return graph(graph_nodes, edges)
    

# define the graph data structure : node, edge and graph

class node():
    def __init__(self,name):
        self.name = name 
        self.edges = set()
        self.partition = None
    # Compute the value of D, D(a) = Ec(a) - Enc(a)
    def compute_D(self):
        Ec = sum(1 for edge in self.edges if edge.get_opposite_node(self).partition != self.partition )
        Enc = len(self.edges) - Ec
        D_ = Ec - Enc
        #print(f"name = {self.name} , Ec = {Ec} , Enc = {Enc} , D_ = {D_}")
        return D_
    
    def add_edge(self,edge):
        self.edges.add(edge)

class edge() : 
    def __init__ (self, u,v):
        self.u = u
        self.v = v
    # Find the node on the opposite side of the edge
    def get_opposite_node(self, node):
        if node.name == self.u :
            return node.graph.get_node(self.v)
        elif node.name == self.v :
            return node.graph.get_node(self.u)

class graph:
    def __init__(self, nodes, edges):
        self.nodes = nodes
        self.edges = edges

        for node in self.nodes:
            node.graph = self

        for edge in self.edges:
            U_ = self.get_node(edge.u)
            V_ = self.get_node(edge.v)
            U_.add_edge(edge)
            V_.add_edge(edge)
    
    # Used to compute the total cost of the graph
    def compute_total_cost(self):
        return sum(1 for edge in self.edges if edge.get_opposite_node(self.get_node(edge.u)).partition != edge.get_opposite_node(self.get_node(edge.v)).partition)

    # Get the node with given node_id
    def get_node(self, node_id):
        return next((v for v in self.nodes if v.name == node_id), None)
    
    # Get the list of nodes in partition
    def get_nodes_in_partition(self, partition):
        return [v for v in self.nodes if v.partition == partition]
    
    # printing the cut set
    def print_cut_set(self):
        print("Cut set: ", [(edge.u,edge.v) for edge in self.edges if edge.get_opposite_node(self.get_node(edge.u)).partition != edge.get_opposite_node(self.get_node(edge.v)).partition])


# Main class to perform the KL bi-parititioning
class KL_bipartition:

    def __init__(self, graph, partition1_name, partition2_name):
        self.graph = graph
        self.partition1_name = partition1_name
        self.partition2_name = partition2_name
        
    # Partition the graph 
    def partition(self):

        # Initialize the partitions - {1 : n/2} and {n/2 + 1 : n}, name them for tracking
        for i, node in enumerate(self.graph.nodes):
            node.partition = self.partition1_name if i < len(self.graph.nodes) // 2 else self.partition2_name
            #print(node.partition)

        p = 0  # pass number
        total_gain = 0

        # Check if the gain can be reduced further, start with the swapped partition
        # Stops when the gain is negative or 0
        while True:
            #debug#print("-------------------------Start of New Pass-------------------------")
            partition0 = [v for v in self.graph.nodes if v.partition == self.partition1_name]
            partition1 = [v for v in self.graph.nodes if v.partition == self.partition2_name]
            #self.graph.print_cut_set()
            #debug#print(f"Partition 0: {[node.name for node in partition0]}\t")
            #debug#print(f"Partition 1: {[node.name for node in partition1]}\t")

            D_ = {v.name: v.compute_D() for v in self.graph.nodes}  # the D value of a node, list comprehension
            G_i = []

            # Run until all the nodes are swapped - final iteration gives 0 gain : stop
            for _ in range(len(self.graph.nodes) // 2):
                max_gain = float('-inf')
                fixed_pair = None

                """for a in partition0 : 
                    print(f"D({a.name}) = {D_[a.name]}")
                for b in partition1 : 
                    print(f"D_({b.name}) = {D_[b.name]}")"""

                # Compute and update the maximum gain obtained in a specific partition configuration
                for a in partition0:
                    for b in partition1:
                        c_ab = len(a.edges & b.edges)
                        gain = D_[a.name] + D_[b.name] - 2 * c_ab

                        if gain > max_gain:
                            max_gain = gain
                            fixed_pair = (a, b)
                            #print(f"gain obtained {max_gain} and fixed_pair = {a.name},{b.name}")

                if not fixed_pair:
                    #print("No selected pair, breaking")
                    break
                
                #print(f"Before update: gain obtained {max_gain} and selected pair = {fixed_pair[0].name},{fixed_pair[1].name}")
                # Fix the nodes, and remove from the working partition  i.e. which nodes can be swapped
                a, b = fixed_pair
                partition0.remove(a)
                partition1.remove(b)
                # Add the nodes to a list, to be used later
                G_i.append((a, b, max_gain))

                # Compute the new D_ values after fixing the pair
                for x in partition0:
                    c_xa = len(x.edges & a.edges)
                    c_xb = len(x.edges & b.edges)
                    D_[x.name] += 2 * c_xa - 2 * c_xb

                for y in partition1:
                    c_yb = len(y.edges & b.edges)
                    c_ya = len(y.edges & a.edges)
                    D_[y.name] += 2 * c_yb - 2 * c_ya

            # Stepwise gains have been calculated : Now calculate the transition that gives max cumulative gain
            G_max = float('-inf')
            j_max = 0
            cumulative_gain = 0

            for j, (a, b, gain) in enumerate(G_i):
                cumulative_gain += gain
                if cumulative_gain > G_max:
                    G_max = cumulative_gain
                    j_max = j + 1
            #print(f"gain = {G_max}")
            if G_max > 0:
                for a, b, _ in G_i[:j_max]:
                    a.partition, b.partition = b.partition, a.partition
                    #debug#print(f"Pass {p} |\tSwapped : ({a.name},{b.name})")
                total_gain += G_max
                #print(f"Pass {p} |\tGain: {G_max}")
                p += 1
                #print("-------------------------Moving To Next Pass-------------------------")
            else:
                break

# Main function

def main(input_filename):
    graph = read_file(input_filename)
    start_time = time.time()
    kl = KL_bipartition(graph,"P0","P1")
    kl.partition()
    end_time = time.time()
    execution_time = end_time - start_time
    #print("---------------------------Run Complete---------------------------")
    #print("-------------------------Printing Summary-------------------------")
    #print("Execution time: {:.2f} seconds".format(execution_time))
    print("Final Partition 0:", [v.name for v in graph.nodes if v.partition == "P0"])
    print("Final Partition 1:", [v.name for v in graph.nodes if v.partition == "P1"])
    #print("Cut Set = ", [ (e.u, e.v) for e in graph.nodes if  ])
    graph.print_cut_set()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 das2541.py input_filename")
        sys.exit(1)
    input_filename = sys.argv[1]
    main(input_filename)