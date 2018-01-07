import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.List;

/**
 * Created by Pati on 16.12.2017.
 */
public class MedianClusterProblem50 {

    CreateGraph50 createGraph = new CreateGraph50("/Users/Pati/Downloads/iris.dat");

    private double[][] initialGraph;
    private double[][] mst;
    private double[][] treeWithoutBiggestValues;

    private List<Integer> seperatedVertices = new ArrayList<>();
    private Set<Integer> cluster1;
    private Set<Integer> cluster2;
    private Set<Integer> cluster3;
    private double sumOfWeights;

    public static void main(String[] args) {

        MedianClusterProblem50 mc = new MedianClusterProblem50();
        mc.solveMedianClusterProblem();

    }

    private void solveMedianClusterProblem() {
        int k = 3; //number of clusters

        //Calculating MST from the initial graph
        initialGraph = createGraph.getDMatrix();
        mst = this.PrimMST(initialGraph);

        //Removing k-1 edges with the biggest weights
        treeWithoutBiggestValues = this.findKMaxAndRemove(mst, k);
//        System.out.println(seperatedVertices);

        int m = treeWithoutBiggestValues.length;

        //calculating k clusters with dfs
        //this code is strictly for k=3

        boolean visited[] = new boolean[m];
        this.dfs(treeWithoutBiggestValues, seperatedVertices.get(0), visited);
        cluster1 = new HashSet<>();

        for (int i = 0; i < visited.length; i++) {
            if (visited[i] == true) {
                cluster1.add(i);
            }
        }
        System.out.println("Cluster 1: " + cluster1);

        this.dfs(treeWithoutBiggestValues, seperatedVertices.get(1), visited);
        cluster2 = new HashSet<>();

        for (int i = 0; i < m; i++) {
            if (visited[i] == true) {
                if (!cluster1.contains(i)) cluster2.add(i);
            }
        }
        System.out.println("Cluster 2: " + cluster2);

        cluster3 = new HashSet<>();
        for (int i = 0; i < m; i++) {
            if (!visited[i]) cluster3.add(i);
        }

        System.out.println("Cluster 3: " + cluster3);

        //Calculating centroids for each cluster and sum of distances to it
        List<Integer> cluster1vertices = new ArrayList<>();
        List<Double> cluster1distances = new ArrayList<>();
        for (int i : cluster1) {
            double distance = 0;
            cluster1vertices.add(i);
            for (int j : cluster1) {
                distance += initialGraph[i][j];
            }
            cluster1distances.add(distance);
        }

        double min1 = Double.POSITIVE_INFINITY;
        int CentrVertex1 = 0;
        for (int i = 0; i < cluster1distances.size(); i++) {
            double d = cluster1distances.get(i);
            if (d < min1) {
                min1 = d;
                CentrVertex1 = cluster1vertices.get(i);
            }
        }
        System.out.println("Centroid of cluster1: " + CentrVertex1 +
                " with a sum of distances to it: " + min1);

        List<Integer> cluster2vertices = new ArrayList<>();
        List<Double> cluster2distances = new ArrayList<>();
        for (int i : cluster2) {
            double distance = 0;
            cluster2vertices.add(i);
            for (int j : cluster2) {
                distance += initialGraph[i][j];
            }
            cluster2distances.add(distance);
        }

        double min2 = Double.POSITIVE_INFINITY;
        int CentrVertex2 = 0;
        for (int i = 0; i < cluster2distances.size(); i++) {
            double d = cluster2distances.get(i);
            if (d < min2) {
                min2 = d;
                CentrVertex2 = cluster2vertices.get(i);
            }
        }
        System.out.println("Centroid of cluster2: " + CentrVertex2 +
                " with a sum of distances to it: " + min2);


        //objective function
        System.out.println("Objective value: " + (min1 + min2));

    }

    //Prim's algorithm
    private double[][] PrimMST(double[][] initialGraph) {
        int m = initialGraph.length;
        sumOfWeights = 0;
        final Set<Integer> unvisited = new HashSet<>();
        final Set<Integer> visited = new HashSet<>();

        for (int i = 0; i < m; i++) {
            unvisited.add(i);
        }
        int startVertex = 0;
        unvisited.remove(startVertex);
        visited.add(startVertex);

        mst = new double[m][m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                mst[i][j] = 0;
            }
        }

        int vertexTo = 0;
        while (!unvisited.isEmpty()) {
            double min = Double.POSITIVE_INFINITY;
            for (int i : visited) {
                for (int j = 0; j < m; j++) {
                    if (unvisited.contains(j)) {
                        double weight = initialGraph[i][j];
                        if (weight < min) {
                            min = weight;
                            startVertex = i;
                            vertexTo = j;
                        }
                    }
                }
            }
            visited.add(vertexTo);
            unvisited.remove(vertexTo);

            mst[startVertex][vertexTo] = mst[vertexTo][startVertex] = min;
            startVertex = vertexTo;
            sumOfWeights += min;
        }
        return mst;
    }

    //Find and remove k-1 edges with the biggest weights
    private double[][] findKMaxAndRemove(double[][] matrix, int k) {
        for (int x = 0; x < k-1; x++) {

            double l = matrix[0][0];
            int vertexFrom = 0;
            int vertexTo = 0;
            for (int i = 0; i < matrix.length; i++) {
                for (int j = 0; j < matrix.length; j++)
                    if (matrix[i][j] > l) {
                        l = matrix[i][j];
                        vertexFrom = i;
                        vertexTo = j;
                    }
            }
            seperatedVertices.add(vertexFrom);
            seperatedVertices.add(vertexTo);

            matrix[vertexFrom][vertexTo] = matrix[vertexTo][vertexFrom] = 0;
            sumOfWeights = sumOfWeights - l;
        }

        return matrix;
    }

    //Depth-First Search for finding paths from vertex s
    void dfs(double[][] matrix, int s, boolean visited[]) {
        visited[s] = true;
        for (int i = 0; i < matrix.length; i++) {
            if (matrix[s][i] != 0 && !visited[i]) {
                dfs(matrix, i, visited);
            }
        }
    }

}

class CreateGraph50 {

    private Scanner sc;

    // not really needed to solve the problem but in case...
    private List<String> pointNames = new ArrayList<>();
    private List<List<Double>> points = new ArrayList<>();

    public double[][] getDMatrix() {
        return DMatrix;
    }

    private double[][] DMatrix;

    int m;
    int n;

    public CreateGraph50(String filePath) {
        try {
            sc = new Scanner(new File(filePath));
        } catch (FileNotFoundException e) {
            System.out.println("Error: Could not find the file specified.");
            System.exit(0);
        }

        initialize();

    }


    public void initialize() {
        m = 50;
        n = 4;

        for (int i = 0; i < m; i++) {

            String pointName = sc.next().substring(3);
            pointNames.add(pointName);
            List<Double> pointVars = new ArrayList<>();
            for (int j = 0; j < n; j++) {

                Double pointVar = sc.nextDouble();
                pointVars.add(pointVar);
            }
            points.add(pointVars);
        }
        DMatrix = new double[m][m];

        for (int i = 0; i < m; i++) {
            sc.next();
            for (int j = 0; j < m; j++) {

                DMatrix[i][j] = sc.nextDouble();


            }
        }

    }
}

