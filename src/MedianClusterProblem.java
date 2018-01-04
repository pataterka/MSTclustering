import edu.princeton.cs.algs4.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Created by Pati on 15.12.2017.
 */
public class MedianClusterProblem {

    MakeGraph createGraph = new MakeGraph("/Users/Pati/Downloads/iris2.dat");
    private EdgeWeightedGraph initialGraph;
    private Graph graphWithRemovedEdges;
    private EdgeWeightedGraph auxGraph;
    private List<Integer> cluster1 = new ArrayList<>();
    private List<Integer> cluster2 = new ArrayList<>();
    private List<Integer> cluster3 = new ArrayList<>();
    private double obj;

    private void solveMedianClusterProblem() {
        int k = 3;
        initialGraph = createGraph.getGraph();
        PrimMST prim = new PrimMST(initialGraph);
        auxGraph = new EdgeWeightedGraph(initialGraph.V());

        List<Edge> sortedEdges = new ArrayList<Edge>();


        for (Edge e : prim.edges()) {

//            auxGraph.addEdge(e);
            sortedEdges.add(e);

        }

        Collections.sort(sortedEdges);

        //this piece of code is implemented just for k=3
        // vertices 1 and 2 are the endpoints of deleted edge, hence they belong to different clusters

        int vertex1 = sortedEdges.get(sortedEdges.size() - 1).either();
        cluster1.add(vertex1);
        int vertex2 = sortedEdges.get(sortedEdges.size() - 1).other(vertex1);
        cluster2.add(vertex2);

        for (int i = 0; i < k - 1; i++) {
            sortedEdges.remove(sortedEdges.size() - 1);
        }
        graphWithRemovedEdges = new Graph(initialGraph.V());
        //obj = null;

        for (Edge e : sortedEdges) {
            obj = obj + e.weight();
            int v = e.either();
            int w = e.other(v);

            graphWithRemovedEdges.addEdge(v, w);
        }

        DepthFirstPaths dfs1 = new DepthFirstPaths(graphWithRemovedEdges, vertex1);

        for (int i = 0; i < initialGraph.V(); i++) {
            if (i == vertex1) continue;
            else if (dfs1.hasPathTo(i)) cluster1.add(i);
        }

        DepthFirstPaths dfs2 = new DepthFirstPaths(graphWithRemovedEdges, vertex2);

        for (int i = 0; i < initialGraph.V(); i++) {
            if (i == vertex2) continue;
            else if (dfs2.hasPathTo(i)) cluster2.add(i);
        }

        for (int i = 0; i < initialGraph.V(); i++) {
            if (!cluster1.contains(i) && !cluster2.contains(i)) cluster3.add(i);
        }


        System.out.println("Objective value: " + obj);
        System.out.println(cluster1);
        System.out.println(cluster2);
        System.out.println(cluster3);

//        List<Integer> cluster1vertices = new ArrayList<>();
//        List<Double> cluster1distances = new ArrayList<>();
//        for (int i : cluster1) {
//            double distance = 0;
//            cluster1vertices.add(i);
//            for (int j : cluster1) {
//                distance += graphWithRemovedEdges.;
//            }
//            cluster1distances.add(distance);
//        }
//        double min1 = Double.POSITIVE_INFINITY;
//        int CentrVertex1 = 0;
//        for (int i = 0; i < cluster1distances.size(); i++) {
//            double d = cluster1distances.get(i);
//            if (d < min1) {
//                min1 = d;
//                CentrVertex1 = cluster1vertices.get(i);
//            }
//        }
//        System.out.println("Centroid of cluster1: " + CentrVertex1 +
//                " with a sum of distances to it: " + min1);
//
//        List<Integer> cluster2vertices = new ArrayList<>();
//        List<Double> cluster2distances = new ArrayList<>();
//        for (int i : cluster2) {
//            double distance = 0;
//            cluster2vertices.add(i);
//            for (int j : cluster2) {
//                distance += initialGraph[i][j];
//            }
//            cluster2distances.add(distance);
//        }
//
//        double min2 = Double.POSITIVE_INFINITY;
//        int CentrVertex2 = 0;
//        for (int i = 0; i < cluster2distances.size(); i++) {
//            double d = cluster2distances.get(i);
//            if (d < min2) {
//                min2 = d;
//                CentrVertex2 = cluster2vertices.get(i);
//            }
//        }
//        System.out.println("Centroid of cluster2: " + CentrVertex2 +
//                " with a sum of distances to it: " + min2);
//
//        List<Integer> cluster3vertices = new ArrayList<>();
//        List<Double> cluster3distances = new ArrayList<>();
//        for (int i : cluster2) {
//            double distance = 0;
//            cluster3vertices.add(i);
//            for (int j : cluster3) {
//                distance += initialGraph[i][j];
//            }
//            cluster3distances.add(distance);
//        }
//
//        double min3 = Double.POSITIVE_INFINITY;
//        int CentrVertex3 = 0;
//        for (int i = 0; i < cluster3distances.size(); i++) {
//            double d = cluster3distances.get(i);
//            if (d < min3) {
//                min3 = d;
//                CentrVertex3 = cluster3vertices.get(i);
//            }
//        }
//        System.out.println("Centroid of cluster3: " + CentrVertex3 +
//                " with a sum of distances to it: " + min3);
//
//
//        //objective function
//        System.out.println("Objective value: " + (min1 + min2 + min3));
    }



    public static void main(String[] args) {

        MedianClusterProblem mc = new MedianClusterProblem();
        mc.solveMedianClusterProblem();

    }

}

class MakeGraph {

    private Scanner sc;


    // not really needed but in case...
    private List<String> pointNames = new ArrayList<>();
    private List<List<Double>> points = new ArrayList<>();
    private EdgeWeightedGraph graph;
    public MakeGraph(String filePath) {
        try {
            sc = new Scanner(new File(filePath));
        } catch (FileNotFoundException e) {
            System.out.println("Error: Could not find the file specified.");
            System.exit(0);
        }

        initialize();

    }

    public EdgeWeightedGraph getGraph() {
        return graph;
    }

    public void initialize() {
        int m = 150;
        int n = 4;

        graph = new EdgeWeightedGraph(m);


        for (int i = 0; i< m; i++){
            sc.next();

//            String pointName = sc.next().substring(3);
//            pointNames.add(pointName);
//            List<Double> pointVars = new ArrayList<>();
//            for (int j = 0; j < n; j++){
//
//                Double pointVar = sc.nextDouble();
//                pointVars.add(pointVar);
//            }
//            points.add(pointVars);
        }

        for (int i = 0; i<m; i++){
            sc.next();
            for (int j = 0; j<m; j++){

                if (i==j){
                    sc.nextDouble();
                }
                else{
                    Edge edge = new Edge(i, j, sc.nextDouble());
                    graph.addEdge(edge);}

            }
        }


    }
}
