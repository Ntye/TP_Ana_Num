package LaplaceFiniteDifference2D;

import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;


public class LaplaceSolverGUI {

    // --- Définitions copiées/adaptées de ODEFiniteDifference2D ---
    interface UValueProvider2D {
        double getValue(double x, double y);
        double getLaplacian(double x, double y); // Pour Δu
        String getName();
    }

    static class ExactSolutionLaplace implements UValueProvider2D {
        public double getValue(double x, double y) {
            return 2 * x + y;
        }
        public double getLaplacian(double x, double y) {
            return 0.0; // Δ(2x+y) = 0
        }
        public String getName() { return "u(x,y) = 2x + y"; }
    }

    static class Solution2D {
        double[][] x_coords_node;
        double[][] y_coords_node;
        double[][] numerical;
        double[][] analytical;
        double errorLinf;
        int nx, ny;
        UValueProvider2D exactSolutionProvider;

        Solution2D(int nx_intervals, int ny_intervals, UValueProvider2D exactSolProvider) {
            this.nx = nx_intervals;
            this.ny = ny_intervals;
            this.x_coords_node = new double[nx + 1][ny + 1];
            this.y_coords_node = new double[nx + 1][ny + 1];
            this.numerical = new double[nx + 1][ny + 1];
            this.analytical = new double[nx + 1][ny + 1];
            this.exactSolutionProvider = exactSolProvider;

            double hx_step = 1.0 / nx;
            double hy_step = 1.0 / ny;
            for (int i = 0; i <= nx; i++) {
                for (int j = 0; j <= ny; j++) {
                    x_coords_node[i][j] = i * hx_step;
                    y_coords_node[i][j] = j * hy_step;
                }
            }
        }
    }

    private static void solveJacobi(double[][] u_solution_grid, double[][] f_source_grid,
                                    int nx_intervals, int ny_intervals, int maxIterations,
                                    double convergenceTolerance, UValueProvider2D uExactProvider) {
        double hx_step = 1.0 / nx_intervals;
        double hy_step = 1.0 / ny_intervals;
        double hx_sq = hx_step * hx_step;
        double hy_sq = hy_step * hy_step;

        double[][] u_old_iter = new double[nx_intervals + 1][ny_intervals + 1];
        double maxAbsoluteDifference = 0.0;

        for (int i = 0; i <= nx_intervals; i++) {
            for (int j = 0; j <= ny_intervals; j++) {
                if (i == 0 || i == nx_intervals || j == 0 || j == ny_intervals) {
                    u_solution_grid[i][j] = uExactProvider.getValue(i * hx_step, j * hy_step);
                } else {
                    u_solution_grid[i][j] = 0.0;
                }
            }
        }

        for (int iter = 0; iter < maxIterations; iter++) {
            for (int i = 0; i <= nx_intervals; i++) {
                System.arraycopy(u_solution_grid[i], 0, u_old_iter[i], 0, ny_intervals + 1);
            }
            maxAbsoluteDifference = 0.0;
            for (int i = 1; i < nx_intervals; i++) {
                for (int j = 1; j < ny_intervals; j++) {
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
            if (maxAbsoluteDifference < convergenceTolerance) {
                System.out.println("Jacobi a convergé en " + (iter + 1) + " itérations. Différence Max = " + maxAbsoluteDifference);
                return;
            }
        }
        System.out.println("Jacobi: Nombre max d'itérations atteint sans convergence. Différence Max = " + maxAbsoluteDifference);
    }
    // --- Fin des définitions copiées/adaptées ---

    /**
     * Résout -Δu = 0 pour u_exact(x,y) = 2x + y sur un maillage Nx x Ny.
     */
    public static Solution2D solveLaplace(int N, UValueProvider2D exactSolutionProvider, int maxIter, double tolerance) {
        int nx_intervals = N;
        int ny_intervals = N;

        Solution2D sol = new Solution2D(nx_intervals, ny_intervals, exactSolutionProvider);
        double hx_step = 1.0 / nx_intervals;
        double hy_step = 1.0 / ny_intervals;

        // Pour -Δu = 0, le terme source f(x,y) est 0 partout.
        // Le laplacien de la solution exacte 2x+y est aussi 0.
        // Donc f_source_terms[i][j] = -exactSolutionProvider.getLaplacian(x,y) sera 0.
        double[][] f_source_terms = new double[nx_intervals + 1][ny_intervals + 1]; // Initialisé à 0 par défaut

        for (int i = 0; i <= nx_intervals; i++) {
            for (int j = 0; j <= ny_intervals; j++) {
                double x = i * hx_step;
                double y = j * hy_step;
                sol.analytical[i][j] = exactSolutionProvider.getValue(x,y);
                // f_source_terms[i][j] = 0; // Explicitement, bien que déjà 0 par défaut
            }
        }

        solveJacobi(sol.numerical, f_source_terms, nx_intervals, ny_intervals, maxIter, tolerance, exactSolutionProvider);

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

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            Font uiFont = new Font("SansSerif", Font.PLAIN, 12);
            Font boldFont = new Font("SansSerif", Font.BOLD, 12);
            Font titleFont = new Font("SansSerif", Font.BOLD, 14);

            JFrame mainFrame = new JFrame("Solveur Laplace 2D pour u=2x+y (Maillage Variable)");
            mainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            mainFrame.setLayout(new BorderLayout(10, 10));
            ((JPanel)mainFrame.getContentPane()).setBorder(BorderFactory.createEmptyBorder(10,10,10,10));

            // Panneau de contrôle
            JPanel controlPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
            controlPanel.add(new JLabel("Nombre d'intervalles N:"));
            JSpinner nSpinner = new JSpinner(new SpinnerNumberModel(10, 2, 100, 1)); // N=10 par défaut, min 2, max 100
            nSpinner.setFont(uiFont);
            controlPanel.add(nSpinner);

            JButton solveButton = new JButton("Résoudre");
            solveButton.setFont(uiFont);
            controlPanel.add(solveButton);

            // Zone pour afficher les résultats textuels
            JTextArea resultsArea = new JTextArea(8, 45);
            resultsArea.setFont(new Font("Monospaced", Font.PLAIN, 12));
            resultsArea.setEditable(false);

            // Panneau pour les heatmaps
            JPanel heatmapsOuterPanel = new JPanel(new BorderLayout(5,5));
            JLabel heatmapTitleLabel = new JLabel("Cartes de Chaleur", SwingConstants.CENTER);
            heatmapTitleLabel.setFont(titleFont);
            heatmapsOuterPanel.add(heatmapTitleLabel, BorderLayout.NORTH);
            JPanel heatmapsGridPanel = new JPanel(new GridLayout(1,3,5,5));
            heatmapsOuterPanel.add(heatmapsGridPanel, BorderLayout.CENTER);

            // Organisation générale
            JPanel textAndHeatmapPanel = new JPanel(new BorderLayout(5,10));
            textAndHeatmapPanel.add(new JScrollPane(resultsArea), BorderLayout.NORTH);
            textAndHeatmapPanel.add(heatmapsOuterPanel, BorderLayout.CENTER);

            mainFrame.add(controlPanel, BorderLayout.NORTH);
            mainFrame.add(textAndHeatmapPanel, BorderLayout.CENTER);

            solveButton.addActionListener(e -> {
                int N = (Integer) nSpinner.getValue();
                UValueProvider2D exactLaplaceSolution = new ExactSolutionLaplace();

                Solution2D sol = solveLaplace(N, exactLaplaceSolution, 20000, 1e-9);

                resultsArea.setText(""); // Effacer les résultats précédents
                StringBuilder sb = new StringBuilder();
                sb.append(String.format("Équation: -Δu = 0 (Solution exacte u=2x+y)\n"));
                sb.append(String.format("Maillage: %dx%d intervalles (h=%.3f)\n", sol.nx, sol.ny, 1.0/sol.nx));
                sb.append(String.format("Erreur L∞: %.4e\n", sol.errorLinf));

                // Si N est petit, afficher la grille dans la zone de texte
                if (N <= 5) { // Par exemple, pour N petit
                    sb.append("\nGrille de solution numérique u_ij:\n");
                    for (int j = sol.ny; j >= 0; j--) {
                        for (int i = 0; i <= sol.nx; i++) {
                            sb.append(String.format("%.3f\t", sol.numerical[i][j]));
                        }
                        sb.append("\n");
                    }
                }
                resultsArea.setText(sb.toString());

                heatmapsGridPanel.removeAll();
                heatmapsGridPanel.add(new HeatmapPanel(sol.numerical, "Numérique (N=" + sol.nx + ")"));
                heatmapsGridPanel.add(new HeatmapPanel(sol.analytical, "Analytique (N=" + sol.nx + ")"));

                double[][] errorGrid = new double[sol.nx+1][sol.ny+1];
                for(int i=0; i<=sol.nx; i++) for(int j=0; j<=sol.ny; j++) errorGrid[i][j] = Math.abs(sol.numerical[i][j] - sol.analytical[i][j]);
                heatmapsGridPanel.add(new HeatmapPanel(errorGrid, "Erreur Absolue (N=" + sol.nx + ")"));

                heatmapsGridPanel.revalidate();
                heatmapsGridPanel.repaint();
                mainFrame.pack();
            });

            mainFrame.pack();
            mainFrame.setMinimumSize(new Dimension(600,500));
            mainFrame.setLocationRelativeTo(null);
            mainFrame.setVisible(true);

            // Lancer une résolution initiale
            solveButton.doClick();
        });
    }
}

// Classe HeatmapPanel (copiée de ODEFiniteDifference2.java et traduite)
class HeatmapPanel extends JPanel {
    private double[][] data;
    private String title;
    private double minVal = Double.MAX_VALUE, maxVal = Double.MIN_VALUE;

    public HeatmapPanel(double[][] dataGrid, String panelTitle) {
        this.data = dataGrid;
        this.title = panelTitle;
        setPreferredSize(new Dimension(220, 250));

        if (data == null || data.length == 0 || data[0].length == 0) return;

        for (double[] row : data) {
            for (double val : row) {
                if (val < minVal) minVal = val;
                if (val > maxVal) maxVal = val;
            }
        }
        if (Math.abs(maxVal - minVal) < 1e-9) {
            maxVal = minVal + 0.5;
            minVal = minVal - 0.5;
        }
        if (minVal == maxVal) maxVal = minVal + 1e-9;
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        if (data == null || data.length == 0 || data[0].length == 0) {
            g.setFont(new Font("SansSerif", Font.PLAIN, 12));
            g.drawString("Aucune donnée à afficher", 10, 20);
            return;
        }

        int panelWidth = getWidth();
        int panelHeight = getHeight();
        int usableWidth = panelWidth - 20;
        int usableHeight = panelHeight - 40;

        int rows = data.length;
        int cols = data[0].length;

        int cellWidth = Math.max(1, usableWidth / cols);
        int cellHeight = Math.max(1, usableHeight / rows);

        int offsetX = (panelWidth - cols * cellWidth) / 2;
        int offsetY = 20 + (usableHeight - rows * cellHeight) / 2;

        g.setColor(Color.BLACK);
        g.setFont(new Font("SansSerif", Font.BOLD, 12));
        g.drawString(title, panelWidth/2 - g.getFontMetrics().stringWidth(title)/2, 15);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                double value = data[i][j];
                float normalized = 0.5f;
                if (maxVal > minVal) {
                     normalized = (float) ((value - minVal) / (maxVal - minVal));
                }
                normalized = Math.max(0f, Math.min(1f, normalized));

                Color color;
                if (normalized < 0.5f) {
                    color = new Color(normalized * 2, normalized * 2, 1f);
                } else {
                     color = new Color(1f, (1 - normalized) * 2, (1 - normalized) * 2);
                }
                g.setColor(color);
                g.fillRect(offsetX + j * cellWidth, offsetY + (rows - 1 - i) * cellHeight, cellWidth, cellHeight);
            }
        }
        g.setColor(Color.BLACK);
        g.setFont(new Font("SansSerif", Font.PLAIN, 10));
        g.drawString(String.format("Min: %.2e", minVal), offsetX, panelHeight - 5);
        g.drawString(String.format("Max: %.2e", maxVal), offsetX + cols * cellWidth - g.getFontMetrics().stringWidth(String.format("Max: %.2e", maxVal)), panelHeight - 5);
    }
}
