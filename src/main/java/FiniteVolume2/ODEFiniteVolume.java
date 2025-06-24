package FiniteVolume2;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Path2D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class ODEFiniteVolume {

    // Interface pour la solution exacte u(x,y) et son laplacien
    interface UValueProvider2D {
        double getValue(double x, double y);
        double getLaplacian(double x, double y);
        String getName();
    }

    // Implémentations de UValueProvider2D
    static class USinPiXSinPiY implements UValueProvider2D {
        public double getValue(double x, double y) {
            return Math.sin(Math.PI * x) * Math.sin(Math.PI * y);
        }
        public double getLaplacian(double x, double y) {
            return -2 * Math.PI * Math.PI * Math.sin(Math.PI * x) * Math.sin(Math.PI * y);
        }
        public String getName() { return "u(x,y) = sin(πx)sin(πy)"; }
    }

    static class UX3Y3 implements UValueProvider2D {
        public double getValue(double x, double y) {
            return x * x * x * y * y * y;
        }
        public double getLaplacian(double x, double y) {
            return 6 * x * y * y * y + 6 * x * x * x * y;
        }
        public String getName() { return "u(x,y) = x³y³"; }
    }

    // Classe pour stocker les résultats de la solution 2D
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

    // Solveur itératif de Jacobi pour le système 2D résultant de -Δu = f
    private static void solveJacobi(double[][] u_solution_grid, double[][] f_source_grid,
                                    int nx_intervals, int ny_intervals, int maxIterations,
                                    double convergenceTolerance, UValueProvider2D uExactProvider) {
        double hx_step = 1.0 / nx_intervals;
        double hy_step = 1.0 / ny_intervals;
        double hx_sq = hx_step * hx_step;
        double hy_sq = hy_step * hy_step;

        double[][] u_old_iter = new double[nx_intervals + 1][ny_intervals + 1];
        double maxAbsoluteDifference = 0.0;

        // Appliquer les conditions aux limites de Dirichlet initiales
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
                System.out.println("Jacobi (Volumes Finis 2D) a convergé en " + (iter + 1) + " itérations. Différence Max = " + maxAbsoluteDifference);
                return;
            }
        }
        System.out.println("Jacobi (Volumes Finis 2D): Nombre max d'itérations atteint sans convergence. Différence Max = " + maxAbsoluteDifference);
    }

    // Méthode principale de résolution pour le problème 2D
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
                f_source_terms[i][j] = -uExactProvider.getLaplacian(x, y);
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
        if (error1 <= 1e-14 && error2 <= 1e-14) return Double.NaN;
        if (error2 <= 1e-14) return Double.POSITIVE_INFINITY;
        if (error1 <= 1e-14) return Double.NEGATIVE_INFINITY;

        double ratio_err = error1 / error2;
        double ratio_n = (double)n2_intervals / n1_intervals;
        if (ratio_err <= 0 || ratio_n <=0) return Double.NaN;
        return Math.log(ratio_err) / Math.log(ratio_n);
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            UValueProvider2D[] exactSolutions = {new USinPiXSinPiY(), new UX3Y3()};
            Font uiFont = new Font("SansSerif", Font.PLAIN, 12);

            JFrame mainFrame = new JFrame("Résolveur Volumes Finis 2D (-Δu = f)");
            mainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            mainFrame.setLayout(new BorderLayout(5,5));

            JPanel controlPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 5, 5));
            JComboBox<UValueProvider2D> exactSolutionCombo = new JComboBox<>(exactSolutions);
            exactSolutionCombo.setFont(uiFont);
            exactSolutionCombo.setRenderer(new DefaultListCellRenderer() {
                 @Override
                public Component getListCellRendererComponent(JList<?> list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
                    super.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);
                    if (value instanceof UValueProvider2D) setText(((UValueProvider2D) value).getName());
                    setFont(uiFont);
                    return this;
                }
            });

            JLabel exactSolLabel = new JLabel("Sol. exacte u(x,y):");
            exactSolLabel.setFont(uiFont);
            controlPanel.add(exactSolLabel);
            controlPanel.add(exactSolutionCombo);
            
            JButton solveButton = new JButton("Résoudre et Afficher");
            solveButton.setFont(uiFont);
            controlPanel.add(solveButton);

            JTextArea resultsArea = new JTextArea(12, 50);
            resultsArea.setEditable(false);
            resultsArea.setFont(new Font("Monospaced", Font.PLAIN, 12));
            JScrollPane scrollPane = new JScrollPane(resultsArea);

            JPanel heatmapsOuterPanel = new JPanel(new BorderLayout(5,5));
            JLabel heatmapTitleLabel = new JLabel("Cartes de Chaleur (Numérique, Analytique, Erreur) - Vol. Finis 2D", SwingConstants.CENTER);
            heatmapTitleLabel.setFont(new Font("SansSerif", Font.BOLD, 14));
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

                resultsArea.append("Résolution pour solution exacte: " + selectedExactSolution.getName() + "\n");
                resultsArea.append("Méthode: Volumes Finis 2D\n");
                resultsArea.append("-----------------------------------------------------------\n");
                resultsArea.append(String.format("%-5s | %-12s | %-8s\n", "N", "Erreur L∞", "Ordre"));
                resultsArea.append("-----------------------------------------------------------\n");

                for (int N : N_values) {
                    Solution2D sol = solve(N, N, selectedExactSolution, 20000, 1e-10);
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
                    heatmapsGridPanel.add(new HeatmapPanel(lastSol.numerical, "Numérique VF (N=" + lastSol.nx + ")"));
                    heatmapsGridPanel.add(new HeatmapPanel(lastSol.analytical, "Analytique VF (N=" + lastSol.nx + ")"));
                    
                    double[][] errorGrid = new double[lastSol.nx+1][lastSol.ny+1];
                    for(int i=0; i<=lastSol.nx; i++) for(int j=0; j<=lastSol.ny; j++) errorGrid[i][j] = Math.abs(lastSol.numerical[i][j] - lastSol.analytical[i][j]);
                    heatmapsGridPanel.add(new HeatmapPanel(errorGrid, "Erreur Absolue VF (N=" + lastSol.nx + ")"));
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

// Classe simple pour afficher des données 2D sous forme de carte de chaleur
class HeatmapPanel extends JPanel {
    private double[][] data;
    private String title;
    private double minVal = Double.MAX_VALUE, maxVal = Double.MIN_VALUE;

    public HeatmapPanel(double[][] dataGrid, String panelTitle) {
        this.data = dataGrid;
        this.title = panelTitle;
        setPreferredSize(new Dimension(250, 280));

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
