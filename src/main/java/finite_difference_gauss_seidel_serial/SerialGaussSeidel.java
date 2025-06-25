package finite_difference_gauss_seidel_serial;

public class SerialGaussSeidel {

    private int N = 3; // Grid size (NxN)
    private double[][] u; // Solution array
    private double h = 1.0 / (N - 1); // Step size, assuming domain [0,1]x[0,1]

    public SerialGaussSeidel() {
        u = new double[N][N];
        initializeBoundaryConditions();
    }

    private void initializeBoundaryConditions() {
        // u = 2x + y
        // Assuming x and y are indices for simplicity in a 3x3 grid
        // For a more general case, x = i*h, y = j*h
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

    public void solve(int maxIterations, double tolerance) {
        for (int k = 0; k < maxIterations; k++) {
            double maxDiff = 0.0;
            for (int i = 1; i < N - 1; i++) { // Interior points
                for (int j = 1; j < N - 1; j++) { // Interior points
                    double old_u_ij = u[i][j];
                    // For -laplacian u = 0, the update rule is:
                    // u_ij = (u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1}) / 4
                    u[i][j] = (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1]) / 4.0;
                    maxDiff = Math.max(maxDiff, Math.abs(u[i][j] - old_u_ij));
                }
            }
            if (maxDiff < tolerance) {
                System.out.println("Converged after " + (k + 1) + " iterations.");
                return;
            }
        }
        System.out.println("Reached maximum iterations (" + maxIterations + ").");
    }

    public void printSolution() {
        System.out.println("Solution u:");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                System.out.printf("%.4f\t", u[i][j]);
            }
            System.out.println();
        }
    }

    public double[][] getSolution() {
        return u;
    }

    public static void main(String[] args) {
        SerialGaussSeidel solver = new SerialGaussSeidel();
        System.out.println("Initial grid with boundary conditions:");
        solver.printSolution();

        solver.solve(100, 1e-4); // Max 100 iterations, tolerance 1e-4

        System.out.println("\nFinal solution:");
        solver.printSolution();
    }
}
