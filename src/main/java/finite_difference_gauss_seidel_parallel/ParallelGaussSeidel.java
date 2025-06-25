package finite_difference_gauss_seidel_parallel;

import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;
import java.util.concurrent.atomic.DoubleAccumulator;

public class ParallelGaussSeidel {

    private int N = 100; // Grid size (NxN)
    private double[][] u; // Solution array
    private double h = 1.0 / (N - 1); // Step size, assuming domain [0,1]x[0,1]
    private ForkJoinPool forkJoinPool;

    public ParallelGaussSeidel() {
        u = new double[N][N];
        initializeBoundaryConditions();
        // Use a ForkJoinPool with a number of threads equal to available processors
        forkJoinPool = new ForkJoinPool();
    }

    private void initializeBoundaryConditions() {
        // u = 2x + y
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == 0 || i == N - 1 || j == 0 || j == N - 1) { // Boundary points
                    double x = i * h;
                    double y = j * h;
                    u[i][j] = 2 * x + y;
                } else {
                    u[i][j] = 0; // Initial guess for interior points
                }
            }
        }
    }

    // Inner class for parallel computation task
    private class GaussSeidelTask extends RecursiveAction {
        private final int startRow, endRow, color; // 0 for red, 1 for black
        private final DoubleAccumulator maxDiffAccumulator;

        GaussSeidelTask(int startRow, int endRow, int color, DoubleAccumulator maxDiffAccumulator) {
            this.startRow = startRow;
            this.endRow = endRow;
            this.color = color;
            this.maxDiffAccumulator = maxDiffAccumulator;
        }

        @Override
        protected void compute() {
            // Threshold for sequential processing
            if (endRow - startRow < 10) { // Adjust threshold as needed
                for (int i = startRow; i < endRow; i++) {
                    for (int j = 1; j < N - 1; j++) {
                        if ((i + j) % 2 == color) { // Process only cells of the current color
                            double old_u_ij = u[i][j];
                            u[i][j] = (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1]) / 4.0;
                            maxDiffAccumulator.accumulate(Math.abs(u[i][j] - old_u_ij));
                        }
                    }
                }
            } else {
                int midRow = (startRow + endRow) / 2;
                invokeAll(new GaussSeidelTask(startRow, midRow, color, maxDiffAccumulator),
                          new GaussSeidelTask(midRow, endRow, color, maxDiffAccumulator));
            }
        }
    }

    public void solve(int maxIterations, double tolerance) {
        for (int k = 0; k < maxIterations; k++) {
            DoubleAccumulator maxDiffAccumulator = new DoubleAccumulator(Math::max, 0.0);

            // Red pass
            GaussSeidelTask redTask = new GaussSeidelTask(1, N - 1, 0, maxDiffAccumulator);
            forkJoinPool.invoke(redTask);

            // Black pass
            GaussSeidelTask blackTask = new GaussSeidelTask(1, N - 1, 1, maxDiffAccumulator);
            forkJoinPool.invoke(blackTask);

            double maxDiff = maxDiffAccumulator.get();

            if (k % 10 == 0 || maxDiff < tolerance) { // Print progress less frequently for larger grids
                 System.out.printf("Iteration %d, Max Diff: %.6e\n", k + 1, maxDiff);
            }

            if (maxDiff < tolerance) {
                System.out.println("Converged after " + (k + 1) + " iterations.");
                forkJoinPool.shutdown();
                return;
            }
        }
        System.out.println("Reached maximum iterations (" + maxIterations + ").");
        forkJoinPool.shutdown();
    }

    public double[][] getSolution() {
        return u;
    }

    public void printSolution() {
        System.out.println("Solution u (showing corners and center for brevity):");
        System.out.printf("u[0][0] = %.4f\n", u[0][0]);
        System.out.printf("u[0][%d] = %.4f\n", N-1, u[0][N-1]);
        System.out.printf("u[%d][0] = %.4f\n", N-1, u[N-1][0]);
        System.out.printf("u[%d][%d] = %.4f\n", N-1, N-1, u[N-1][N-1]);
        if (N > 2) {
             System.out.printf("u[%d][%d] (center) = %.4f\n", N/2, N/2, u[N/2][N/2]);
        }
    }

    public static void main(String[] args) {
        ParallelGaussSeidel solver = new ParallelGaussSeidel();
        // System.out.println("Initial grid with boundary conditions (corners):");
        // solver.printSolution(); // Printing the whole 100x100 grid is too much

        long startTime = System.currentTimeMillis();
        solver.solve(20000, 1e-6); // Max 20000 iterations, tolerance 1e-6
        long endTime = System.currentTimeMillis();

        System.out.println("\nFinal solution (corners and center):");
        solver.printSolution();
        System.out.println("Time taken: " + (endTime - startTime) + " ms");
    }
}
