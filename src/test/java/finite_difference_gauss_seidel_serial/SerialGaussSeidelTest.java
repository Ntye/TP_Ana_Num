package finite_difference_gauss_seidel_serial;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;

public class SerialGaussSeidelTest {

    @Test
    void testSolve3x3() {
        SerialGaussSeidel solver = new SerialGaussSeidel();
        solver.solve(1000, 1e-6); // Increased iterations for potentially better convergence

        // Accessing u directly for testing is not ideal,
        // but for this specific problem, we'll make an exception or add a getter.
        // For now, let's assume we add a getter getSolution() in SerialGaussSeidel
        // Or, make u package-private or protected for testing.
        // Let's modify SerialGaussSeidel to have a getter for u.

        // Expected value for the center point u[1][1]
        // h = 0.5
        // u[0][1] = 2*(0*0.5) + (1*0.5) = 0.5
        // u[1][0] = 2*(1*0.5) + (0*0.5) = 1.0
        // u[1][2] = 2*(1*0.5) + (2*0.5) = 2.0
        // u[2][1] = 2*(2*0.5) + (1*0.5) = 2.5
        // u[1][1] = (u[0][1] + u[2][1] + u[1][0] + u[1][2]) / 4
        // u[1][1] = (0.5 + 2.5 + 1.0 + 2.0) / 4 = 6.0 / 4 = 1.5

        double[][] solution = solver.getSolution(); // Assuming getter is added
        assertEquals(1.5, solution[1][1], 1e-5);

        // Check boundary values as well
        assertEquals(0.0, solution[0][0], 1e-5);
        assertEquals(0.5, solution[0][1], 1e-5);
        assertEquals(1.0, solution[0][2], 1e-5);
        assertEquals(1.0, solution[1][0], 1e-5);
        assertEquals(2.0, solution[1][2], 1e-5);
        assertEquals(2.0, solution[2][0], 1e-5);
        assertEquals(2.5, solution[2][1], 1e-5);
        assertEquals(3.0, solution[2][2], 1e-5);
    }
}
