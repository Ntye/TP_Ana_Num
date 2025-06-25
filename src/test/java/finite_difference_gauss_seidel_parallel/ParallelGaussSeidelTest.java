package finite_difference_gauss_seidel_parallel;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class ParallelGaussSeidelTest {

    @Test
    void testSolve100x100() {
        ParallelGaussSeidel solver = new ParallelGaussSeidel();
        // Use a sufficiently small tolerance and enough iterations
        // The default in main is 20000 iterations and 1e-6 tolerance.
        // For a test, we might want to be a bit more relaxed or ensure it runs faster.
        // However, the solution is known (u=2x+y), so it should converge well.
        solver.solve(50000, 1e-7); // Increased iterations for test precision

        double[][] solution = solver.getSolution(); // Assuming getter is added
        int N = 100;
        double h = 1.0 / (N - 1);

        // Check boundary values first
        // u[0][0] = 2*(0*h) + (0*h) = 0
        assertEquals(0.0, solution[0][0], 1e-5);
        // u[0][N-1] = 2*(0*h) + ((N-1)*h) = 2*0 + 1 = 1.0
        assertEquals(1.0, solution[0][N-1], 1e-5);
        // u[N-1][0] = 2*((N-1)*h) + (0*h) = 2*1 + 0 = 2.0
        assertEquals(2.0, solution[N-1][0], 1e-5);
        // u[N-1][N-1] = 2*((N-1)*h) + ((N-1)*h) = 2*1 + 1 = 3.0
        assertEquals(3.0, solution[N-1][N-1], 1e-5);

        // Check a few interior points
        // For -laplacian u = 0, and boundary u = 2x + y, the solution is u = 2x + y everywhere.
        int i, j;

        // Point (1,1) (actual indices)
        i = 1; j = 1;
        double x_11 = i * h;
        double y_11 = j * h;
        double expected_u_11 = 2 * x_11 + y_11;
        assertEquals(expected_u_11, solution[i][j], 1e-4);

        // Point (N/4, N/4)
        i = N / 4; j = N / 4;
        double x_N4 = i * h;
        double y_N4 = j * h;
        double expected_u_N4 = 2 * x_N4 + y_N4;
        assertEquals(expected_u_N4, solution[i][j], 1e-4);

        // Point (N/2, N/2) - center
        i = N / 2; j = N / 2;
        double x_N2 = i * h;
        double y_N2 = j * h;
        double expected_u_N2 = 2 * x_N2 + y_N2;
        assertEquals(expected_u_N2, solution[i][j], 1e-4);

        // Point (N-2, N-2) - near the other corner
        i = N - 2; j = N - 2;
        double x_N_2 = i * h;
        double y_N_2 = j * h;
        double expected_u_N_2 = 2 * x_N_2 + y_N_2;
        assertEquals(expected_u_N_2, solution[i][j], 1e-4);

        // Point (20, 70) - arbitrary interior
        i = 20; j = 70;
        double x_20_70 = i*h;
        double y_20_70 = j*h;
        double expected_u_20_70 = 2 * x_20_70 + y_20_70;
        assertEquals(expected_u_20_70, solution[i][j], 1e-4);
    }
}
