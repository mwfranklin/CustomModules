#########################################################################################################################
############################################### GRAPH-ylococcus aureus ##################################################
#########################################################################################################################

class Graph(object):
    ####################################################################
    #### basic graph functions
    ####################################################################
    
    def __init__(self, graph_dict=None):
        """ initializes a graph object 
            If no dictionary or None is given, 
            an empty dictionary will be used
        """
        if graph_dict == None:
            graph_dict = {}
        self.__graph_dict = graph_dict
    
    
    def vertices(self):
        """ returns the vertices of a graph """
        return list(self.__graph_dict.keys())

    
    def edges(self):
        """ returns the edges of a graph """
        return self.__generate_edges()
    
    
    def add_vertex(self, vertex):
        """ If the vertex "vertex" is not in 
            self.__graph_dict, a key "vertex" with an empty
            list as a value is added to the dictionary. 
            Otherwise nothing has to be done. 
        """
        if vertex not in self.__graph_dict:
            self.__graph_dict[vertex] = []

    
    def add_edge(self, vertex1, vertex2):
        """ assumes that edge is of type set, tuple or list; 
            undirected graph
        """
        if vertex1 in self.__graph_dict:
            if vertex2 not in self.__graph_dict[vertex1]:
                self.__graph_dict[vertex1].append(vertex2)
        else:
            self.__graph_dict[vertex1] = [vertex2]
            
        if vertex2 in self.__graph_dict:
            if vertex1 not in self.__graph_dict[vertex2]:
                self.__graph_dict[vertex2].append(vertex1)
        else:
            self.__graph_dict[vertex2] = [vertex1]

            
    def __generate_edges(self):
        """ A static method generating the edges of the 
            graph "graph". Edges are represented as sets 
            with one (a loop back to the vertex) or two 
            vertices 
        """
        edges = []
        for vertex in self.__graph_dict:
            for neighbour in self.__graph_dict[vertex]:
                if {neighbour, vertex} not in edges:
                    edges.append({vertex, neighbour})
        return edges
    
    
        
    ####################################################################
    #### print functions
    ####################################################################
    
    def print_vertices(self):
        vertices_list= list(self.vertices())
        res = "\n\nvertices: "
        for k in vertices_list:
            res += str(k) + "\t"
            if vertices_list.index(k)%10==0: res += "\n"
        print(res)
        
    def print_edges(self):
        edge_list = list(self.__generate_edges())
        res = "\n\nedges:"
        for edge in edge_list:
            res += str(edge) + "\t"
            if edge_list.index(edge)%4==0: res += "\n"
        print(res)
        
        
        
    def print_adjancencey_list(self):
        vertices_list= list(self.vertices())
        res = ''
        for v in vertices_list:
            res += "\nAdjacency list of vertex"
            res += "\nhead "
            for i in self.__graph_dict[v]:
                res += (" -> %s " % i)
        print(res)
    
    ####################################################################
    #### other functions
    ####################################################################
    # helper function for getting list of connected vertices with DFS
    def DFSUtil_connected(self,v,visited):
        # add v to list and mark as visited
        L = [v]
        visited[v] = True
        # Recursive for all the vertices adjacent to this vertex
        for i in self.__graph_dict[v]:
            if visited[i] == False:
                L += self.DFSUtil_connected(i, visited)
        return(L)
    
    ## find a list of all vertices connected to v
    def DFS_connected_graph(self, v):
        # Mark all the vertices as not visited
        visited = dict.fromkeys(self.vertices(),False)
        L = self.DFSUtil_connected(v,visited)
        return(L)
    
    
    # helper function for getting list of connected vertices with BFS
    def BFSUtil_connected(self,v,visited):
        L = []
        S = [v]
        while len(S)>0:
            i = S[0]
            if visited[i] == False:
                L += [i]
                visited[i] = True
                neighbors = list(self.__graph_dict[i])
                for neighbor in neighbors:
                    if visited[neighbor] == False:
                        S += [neighbor]
            S.remove(i)
        return(L)
    
    
    def BFS_connected_graph(self, v):
        # Mark all the vertices as not visited
        visited = dict.fromkeys(self.vertices(),False)
        L = self.BFSUtil_connected(v,visited)
        return(L)
    
    def get_connected_graphs(self):
        graphs = {}
        unchecked_v = list(self.__graph_dict)
        
        graph_num = 0
        while len(unchecked_v)>0:
            graph_num +=1
            graphs[graph_num] = self.BFS_connected_graph(unchecked_v[0])
            
            for v in graphs[graph_num]:
                unchecked_v.remove(v)
        return(graphs)
    
    
    def vertex_degree(self, vertex):
        """ The degree of a vertex is the number of edges connecting
            it, i.e. the number of adjacent vertices. Loops are counted 
            double, i.e. every occurence of vertex in the list 
            of adjacent vertices. """ 
        adj_vertices =  self.__graph_dict[vertex]
        degree = len(adj_vertices) + adj_vertices.count(vertex)
        return degree
    
    
    def degree_sequence(self):
        """ calculates the degree sequence """
        seq = []
        for vertex in self.__graph_dict:
            seq.append(self.vertex_degree(vertex))
        seq.sort(reverse=True)
        return tuple(seq)
    
    def degree_dict(self):
        degree_dict = {}
        """ calculates the degree sequence """
        degree_dict = {}
        for vertex in self.__graph_dict:
            degree_dict[vertex] = self.vertex_degree(vertex)
            
        return degree_dict
    
    
    def find_isolated_vertices(self):
        """ returns a list of isolated vertices. """
        graph = self.__graph_dict
        isolated = []
        for vertex in graph:
            #print(isolated, vertex)
            if not graph[vertex]:
                isolated += [vertex]
        return isolated
    
    
    def neighbors(self, vertex):
        """ returns a list of adjacent vertices. """
        adj_vertices =  self.__graph_dict[vertex]
        return(adj_vertices)