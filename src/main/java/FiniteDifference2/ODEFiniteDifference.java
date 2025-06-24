package FiniteDifference2;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Path2D; // Pourrait être utile pour des tracés plus complexes
import java.awt.image.BufferedImage; // Optionnel pour affichage heatmap
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors; // Ajouté pour Java 8+ streams

public class ODEFiniteDifference {

    // Interface for the exact solution u(x,y) and its Laplacian
    interface UValueProvider2D {
        double getValue(double x, double y); // u(x,y) value
        double getLaplacian(double x, double y); // Laplacian value: u_xx + u_yy
        String getName(); // Name of the exact solution, e.g., "u(x,y) = sin(πx)sin(πy)"
    }

    // Implementations of UValueProvider2D
    static class USinPiXSinPiY implements UValueProvider2D {
        public double getValue(double x, double y) {
            return Math.sin(Math.PI * x) * Math.sin(Math.PI * y);
        }
        public double getLaplacian(double x, double y) {
            // u_xx = -π²sin(πx)sin(πy), u_yy = -π²sin(πx)sin(πy)
            // Δu = -2π²sin(πx)sin(πy)
            return -2 * Math.PI * Math.PI * Math.sin(Math.PI * x) * Math.sin(Math.PI * y);
        }
        public String getName() { return "u(x,y) = sin(πx)sin(πy)"; } // Name for UI
    }

    static class UX3Y3 implements UValueProvider2D {
        public double getValue(double x, double y) {
            return x * x * x * y * y * y;
        }
        public double getLaplacian(double x, double y) {
            // u_x = 3x²y³, u_xx = 6xy³
            // u_y = 3x³y², u_yy = 6x³y
            // Δu = 6xy³ + 6x³y
            return 6 * x * y * y * y + 6 * x * x * x * y;
        }
        public String getName() { return "u(x,y) = x³y³"; } // Name for UI
    }

    // Class to store 2D solution results
    static class Solution2D {
        double[][] x_coords_node; // x_coords_node[i][j] = x_i (x-coordinate of node i,j)
        double[][] y_coords_node; // y_coords_node[i][j] = y_j (y-coordinate of node i,j)
        double[][] numerical;     // Numerical solution u_num[i][j] at node (i,j)
        double[][] analytical;    // Analytical solution u_exact[i][j] at node (i,j)
        double errorLinf;         // L-infinity error
        int nx, ny;               // Number of intervals in x and y directions
        UValueProvider2D exactSolutionProvider; // Provider for the exact solution details

        Solution2D(int nx_intervals, int ny_intervals, UValueProvider2D exactSolProvider) {
            this.nx = nx_intervals;
            this.ny = ny_intervals;
            // Grids have (nx+1) x (ny+1) points/nodes
            this.x_coords_node = new double[nx + 1][ny + 1];
            this.y_coords_node = new double[nx + 1][ny + 1];
            this.numerical = new double[nx + 1][ny + 1];
            this.analytical = new double[nx + 1][ny + 1];
            this.exactSolutionProvider = exactSolProvider;

            double hx_step = 1.0 / nx; // Step size in x
            double hy_step = 1.0 / ny; // Step size in y
            for (int i = 0; i <= nx; i++) {
                for (int j = 0; j <= ny; j++) {
                    x_coords_node[i][j] = i * hx_step;
                    y_coords_node[i][j] = j * hy_step;
                }
            }
        }
    }

    // Jacobi iterative solver for the 2D system from -Δu = f
    // u_solution_grid: 2D grid for the solution (output), includes BCs.
    // f_source_grid: 2D grid for the source term f(x_i, y_j).
    // nx_intervals, ny_intervals: number of intervals.
    private static void solveJacobi(double[][] u_solution_grid, double[][] f_source_grid,
                                    int nx_intervals, int ny_intervals, int maxIterations,
                                    double convergenceTolerance, UValueProvider2D uExactProvider) {
        double hx_step = 1.0 / nx_intervals;
        double hy_step = 1.0 / ny_intervals;
        double hx_sq = hx_step * hx_step;
        double hy_sq = hy_step * hy_step;

        double[][] u_old_iter = new double[nx_intervals + 1][ny_intervals + 1]; // Previous iteration values
        double maxAbsoluteDifference = 0.0; // Declare for scope

        // Apply initial Dirichlet boundary conditions to u_solution_grid
        for (int i = 0; i <= nx_intervals; i++) {
            for (int j = 0; j <= ny_intervals; j++) {
                if (i == 0 || i == nx_intervals || j == 0 || j == ny_intervals) { // On the boundary
                    u_solution_grid[i][j] = uExactProvider.getValue(i * hx_step, j * hy_step);
                } else {
                    u_solution_grid[i][j] = 0.0; // Initial guess for interior points
                }
            }
        }

        for (int iter = 0; iter < maxIterations; iter++) {
            // Copy current solution u_solution_grid to u_old_iter
            for (int i = 0; i <= nx_intervals; i++) {
                System.arraycopy(u_solution_grid[i], 0, u_old_iter[i], 0, ny_intervals + 1);
            }

            maxAbsoluteDifference = 0.0; // Reset max difference for this iteration

            // Update interior points
            for (int i = 1; i < nx_intervals; i++) {
                for (int j = 1; j < ny_intervals; j++) {
                    // Standard 5-point stencil for - (u_xx + u_yy) = f_ij
                    // (2/hx² + 2/hy²) u_ij = (u_old[i+1][j] + u_old[i-1][j])/hx² +
                    //                        (u_old[i][j+1] + u_old[i][j-1])/hy² + f_ij
                    
                    double sum_neighbors_terms = (u_old_iter[i-1][j] + u_old_iter[i+1][j])/hx_sq +
                                                 (u_old_iter[i][j-1] + u_old_iter[i][j+1])/hy_sq;
                    double denominator = (2.0/hx_sq + 2.0/hy_sq);
                    
                    u_solution_grid[i][j] = (sum_neighbors_terms + f_source_grid[i][j]) / denominator;
                    
                    double difference = Math.abs(u_solution_grid[i][j] - u_old_iter[i][j]);
                    if (difference > maxAbsoluteDifference) {
                        maxAbsoluteDifference = difference;
                    }
                }
            }

            // Check for convergence
            if (maxAbsoluteDifference < convergenceTolerance) {
                System.out.println("Jacobi converged in " + (iter + 1) + " iterations. Max Difference = " + maxAbsoluteDifference);
                return;
            }
        }
        System.out.println("Jacobi: Max iterations reached without convergence. Max Difference = " + maxAbsoluteDifference);
    }

    // Main solver method for 2D problem
    public static Solution2D solve(int nx_intervals, int ny_intervals, UValueProvider2D uExactProvider,
                                   int maxSolverIter, double solverTolerance) {
        Solution2D sol = new Solution2D(nx_intervals, ny_intervals, uExactProvider);
        double hx_step = 1.0 / nx_intervals;
        double hy_step = 1.0 / ny_intervals;

        double[][] f_source_terms = new double[nx_intervals + 1][ny_intervals + 1];
        for (int i = 0; i <= nx_intervals; i++) {
            for (int j = 0; j <= ny_intervals; j++) {
                double x = i * hx_step;
                double y = j * hy_step;
                f_source_terms[i][j] = -uExactProvider.getLaplacian(x, y); // f = -Δu_exact
                sol.analytical[i][j] = uExactProvider.getValue(x,y);
            }
        }
        
        solveJacobi(sol.numerical, f_source_terms, nx_intervals, ny_intervals, maxSolverIter, solverTolerance, uExactProvider);

        sol.errorLinf = 0.0;
        for (int i = 0; i <= nx_intervals; i++) {
            for (int j = 0; j <= ny_intervals; j++) {
                double diff = Math.abs(sol.numerical[i][j] - sol.analytical[i][j]);
                if (diff > sol.errorLinf) {
                    sol.errorLinf = diff;
                }
            }
        }
        return sol;
    }
    
    public static double calculateConvergenceOrder(double error1, double error2, int n1_intervals, int n2_intervals) {
        if (n1_intervals == n2_intervals) return Double.NaN;
        if (error1 <= 1e-14 && error2 <= 1e-14) return Double.NaN; // Proche de la précision machine
        if (error2 <= 1e-14) return Double.POSITIVE_INFINITY; // Convergence "parfaite"
        if (error1 <= 1e-14) return Double.NEGATIVE_INFINITY; // Était "parfait", ne l'est plus

        double ratio_err = error1 / error2;
        double ratio_n = (double)n2_intervals / n1_intervals;
        if (ratio_err <= 0 || ratio_n <=0) return Double.NaN; // Arguments de log invalides
        return Math.log(ratio_err) / Math.log(ratio_n);
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            UValueProvider2D[] exactSolutions = {new USinPiXSinPiY(), new UX3Y3()};
            
            JFrame mainFrame = new JFrame("2D Finite Difference Solver (-Δu = f)"); // English title
            mainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            mainFrame.setLayout(new BorderLayout(5,5));

            JPanel controlPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 5, 5));
            JComboBox<UValueProvider2D> exactSolutionCombo = new JComboBox<>(exactSolutions);
            exactSolutionCombo.setRenderer(new DefaultListCellRenderer() {
                 @Override
                public Component getListCellRendererComponent(JList<?> list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
                    super.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);
                    if (value instanceof UValueProvider2D) setText(((UValueProvider2D) value).getName());
                    return this;
                }
            });

            controlPanel.add(new JLabel("Exact Solution u(x,y):")); // English label
            controlPanel.add(exactSolutionCombo);
            
            JButton solveButton = new JButton("Solve & Display"); // English button text
            controlPanel.add(solveButton);

            JTextArea resultsArea = new JTextArea(12, 50);
            resultsArea.setEditable(false);
            resultsArea.setFont(new Font("Monospaced", Font.PLAIN, 12));
            JScrollPane scrollPane = new JScrollPane(resultsArea);

            JPanel heatmapsOuterPanel = new JPanel(new BorderLayout(5,5));
            JLabel heatmapTitleLabel = new JLabel("Heatmaps (Numerical, Analytical, Error) - Finite Differences 2D", SwingConstants.CENTER); // English title
            heatmapsOuterPanel.add(heatmapTitleLabel, BorderLayout.NORTH);
            JPanel heatmapsGridPanel = new JPanel(new GridLayout(1,3,5,5));
            heatmapsOuterPanel.add(heatmapsGridPanel, BorderLayout.CENTER);

            JSplitPane splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, scrollPane, heatmapsOuterPanel);
            splitPane.setResizeWeight(0.4);

            mainFrame.add(controlPanel, BorderLayout.NORTH);
            mainFrame.add(splitPane, BorderLayout.CENTER);

            solveButton.addActionListener(e -> {
                UValueProvider2D selectedExactSolution = (UValueProvider2D) exactSolutionCombo.getSelectedItem();
                resultsArea.setText("");
                heatmapsGridPanel.removeAll();

                int[] N_values = {10, 20, 40, 80};
                List<Solution2D> solutions = new ArrayList<>();
                Solution2D prevSol = null;

                resultsArea.append("Solving for Exact Solution: " + selectedExactSolution.getName() + "\n"); // English
                resultsArea.append("Method: 2D Finite Differences\n");
                resultsArea.append("-----------------------------------------------------------\n");
                resultsArea.append(String.format("%-5s | %-12s | %-8s\n", "N", "L∞ Error", "Order")); // English
                resultsArea.append("-----------------------------------------------------------\n");


                for (int N : N_values) {
                    Solution2D sol = solve(N, N, selectedExactSolution, 20000, 1e-10); // Using high iterations for Jacobi
                    solutions.add(sol);
                    String orderStr = "-";
                    if (prevSol != null) {
                        double order = calculateConvergenceOrder(prevSol.errorLinf, sol.errorLinf, prevSol.nx, sol.nx);
                        orderStr = String.format("%.2f", order);
                    }
                    resultsArea.append(String.format("%-5d | %-12.4e | %-8s\n", N, sol.errorLinf, orderStr));
                    prevSol = sol;
                }
                resultsArea.append("-----------------------------------------------------------\n");

                if (!solutions.isEmpty()) {
                    Solution2D lastSol = solutions.get(solutions.size()-1);
                    heatmapsGridPanel.add(new HeatmapPanel(lastSol.numerical, "Numerical (N=" + lastSol.nx + ")")); // English
                    heatmapsGridPanel.add(new HeatmapPanel(lastSol.analytical, "Analytical (N=" + lastSol.nx + ")"));// English

                    double[][] errorGrid = new double[lastSol.nx+1][lastSol.ny+1];
                    for(int i=0; i<=lastSol.nx; i++) for(int j=0; j<=lastSol.ny; j++) errorGrid[i][j] = Math.abs(lastSol.numerical[i][j] - lastSol.analytical[i][j]);
                    heatmapsGridPanel.add(new HeatmapPanel(errorGrid, "Absolute Error (N=" + lastSol.nx + ")")); // English
                }
                heatmapsGridPanel.revalidate();
                heatmapsGridPanel.repaint();
                mainFrame.pack();
            });
            
            mainFrame.setMinimumSize(new Dimension(600, 700));
            mainFrame.pack();
            mainFrame.setLocationRelativeTo(null);
            mainFrame.setVisible(true);
            
            if (exactSolutionCombo.getItemCount() > 0) {
                 exactSolutionCombo.setSelectedIndex(0);
                 solveButton.doClick();
            }
        });
    }
}

// Simple class to display 2D data as a heatmap
class HeatmapPanel extends JPanel {
    private double[][] data;
    private String title;
    private double minVal = Double.MAX_VALUE, maxVal = Double.MIN_VALUE;

    public HeatmapPanel(double[][] dataGrid, String panelTitle) {
        this.data = dataGrid;
        this.title = panelTitle;
        setPreferredSize(new Dimension(250, 280)); // Preferred size for each heatmap panel

        if (data == null || data.length == 0 || data[0].length == 0) return;

        // Find min and max values in the data for color scaling
        for (double[] row : data) {
            for (double val : row) {
                if (val < minVal) minVal = val;
                if (val > maxVal) maxVal = val;
            }
        }
        if (Math.abs(maxVal - minVal) < 1e-9) { // Handle case where all values are (nearly) equal
            maxVal = minVal + 0.5;
            minVal = minVal - 0.5;
        }
        if (minVal == maxVal) maxVal = minVal + 1e-9; // Ensure maxVal > minVal to avoid division by zero
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        if (data == null || data.length == 0 || data[0].length == 0) {
            g.drawString("No data to display", 10, 20); // English
            return;
        }

        int panelWidth = getWidth();
        int panelHeight = getHeight();
        // Leave space for title and legend
        int usableWidth = panelWidth - 20;
        int usableHeight = panelHeight - 40;

        int rows = data.length;      // Corresponds to ny+1 points (y-direction of grid)
        int cols = data[0].length;   // Corresponds to nx+1 points (x-direction of grid)

        int cellWidth = Math.max(1, usableWidth / cols);
        int cellHeight = Math.max(1, usableHeight / rows);
        
        // Center the heatmap drawing
        int offsetX = (panelWidth - cols * cellWidth) / 2;
        int offsetY = 20 + (usableHeight - rows * cellHeight) / 2; // 20 for title

        g.setColor(Color.BLACK);
        g.drawString(title, panelWidth/2 - g.getFontMetrics().stringWidth(title)/2, 15); // Draw title

        for (int i = 0; i < rows; i++) { // i iterates over rows of `data` (often y-direction)
            for (int j = 0; j < cols; j++) { // j iterates over columns (often x-direction)
                double value = data[i][j];
                float normalized = 0.5f;
                if (maxVal > minVal) {
                     normalized = (float) ((value - minVal) / (maxVal - minVal));
                }
                normalized = Math.max(0f, Math.min(1f, normalized)); // Clamp to [0, 1]

                Color color; // Simple blue (low) to red (high) color scale
                if (normalized < 0.5f) {
                    color = new Color(normalized * 2, normalized * 2, 1f); // Blue -> Cyan -> White-ish
                } else {
                     color = new Color(1f, (1 - normalized) * 2, (1 - normalized) * 2); // White-ish -> Yellow -> Red
                }
                g.setColor(color);
                // Standard display: data[i][j] where i is row index (y), j is col index (x)
                // Graphics origin (0,0) is top-left.
                // To make y-axis point upwards, invert the drawn y-coordinate.
                // (rows - 1 - i) inverts y-axis if data[0] is the bottom row.
                // If data[0] is top row (common in image processing), then just 'i'.
                // Assuming data[i][j] corresponds to u(x_j, y_i) where y_i increases upwards.
                // So, to draw, row i=0 (y=0) should be at the bottom.
                g.fillRect(offsetX + j * cellWidth, offsetY + (rows - 1 - i) * cellHeight, cellWidth, cellHeight);
            }
        }
        g.setColor(Color.BLACK);
        g.drawString(String.format("Min: %.2e", minVal), offsetX, panelHeight - 5);
        g.drawString(String.format("Max: %.2e", maxVal), offsetX + cols * cellWidth - g.getFontMetrics().stringWidth(String.format("Max: %.2e", maxVal)), panelHeight - 5);
    }
}
