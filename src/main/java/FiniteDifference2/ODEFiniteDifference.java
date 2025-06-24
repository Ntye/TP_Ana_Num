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

    // Interface pour la solution exacte u(x,y) et son laplacien
    interface UValueProvider2D {
        double getValue(double x, double y);
        double getLaplacian(double x, double y); // u_xx + u_yy
        String getName();
    }

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

    static class Solution2D {
        double[][] x_coords_node; // x_coords_node[i][j] = x_i (coordonnée x du noeud i,j)
        double[][] y_coords_node; // y_coords_node[i][j] = y_j (coordonnée y du noeud i,j)
        double[][] numerical;     // Valeur numérique u_num[i][j] au noeud (i,j)
        double[][] analytical;    // Valeur analytique u_exact[i][j] au noeud (i,j)
        double errorLinf;
        int nx, ny; // Nombre d'intervalles dans chaque direction (Nx segments en x, Ny en y)
        UValueProvider2D exactSolutionProvider;

        Solution2D(int nx_intervals, int ny_intervals, UValueProvider2D exactSolProvider) {
            this.nx = nx_intervals;
            this.ny = ny_intervals;
            // Les grilles ont (nx+1) x (ny+1) points/noeuds
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
        double maxAbsoluteDifference = 0.0; // Déclarer ici pour la portée

        // Appliquer les conditions aux limites (Dirichlet) à la grille de solution u_solution_grid
        for (int i = 0; i <= nx_intervals; i++) {
            for (int j = 0; j <= ny_intervals; j++) {
                if (i == 0 || i == nx_intervals || j == 0 || j == ny_intervals) { // Sur le bord
                    u_solution_grid[i][j] = uExactProvider.getValue(i * hx_step, j * hy_step);
                } else {
                    u_solution_grid[i][j] = 0.0; // Initialisation pour les points intérieurs (estimation initiale)
                }
            }
        }

        for (int iter = 0; iter < maxIterations; iter++) {
            for (int i = 0; i <= nx_intervals; i++) { // Copier u vers u_old_iter
                System.arraycopy(u_solution_grid[i], 0, u_old_iter[i], 0, ny_intervals + 1);
            }

            maxAbsoluteDifference = 0.0; // Réinitialiser pour cette itération

            // Mettre à jour les points intérieurs
            for (int i = 1; i < nx_intervals; i++) {
                for (int j = 1; j < ny_intervals; j++) {
                    // Pour - ( (u_i+1,j - 2u_ij + u_i-1,j)/hx_sq + (u_i,j+1 - 2u_ij + u_i,j-1)/hy_sq ) = f_ij
                    // (2/hx_sq + 2/hy_sq) u_ij = (u_i+1,j + u_i-1,j)/hx_sq + (u_i,j+1 + u_i,j-1)/hy_sq + f_ij
                    // u_ij = ( (u_old_iter[i+1][j] + u_old_iter[i-1][j])/hx_sq +
                    //          (u_old_iter[i][j+1] + u_old_iter[i][j-1])/hy_sq +
                    //          f_source_grid[i][j]
                    //        ) / (2.0/hx_sq + 2.0/hy_sq)
                    
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
                System.out.println("Jacobi a convergé en " + (iter + 1) + " itérations. MaxDiff = " + maxAbsoluteDifference);
                return;
            }
        }
        System.out.println("Jacobi: convergence non atteinte après " + maxIterations + " itérations. MaxDiff = " + maxAbsoluteDifference);
    }

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
            
            JFrame mainFrame = new JFrame("Résolveur Différences Finies 2D (-Δu = f)");
            mainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            mainFrame.setLayout(new BorderLayout(5,5)); // Ajout de marges

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

            controlPanel.add(new JLabel("Sol. exacte u(x,y):"));
            controlPanel.add(exactSolutionCombo);
            
            JButton solveButton = new JButton("Résoudre et Afficher");
            controlPanel.add(solveButton);

            JTextArea resultsArea = new JTextArea(12, 50); // Taille ajustée
            resultsArea.setEditable(false);
            resultsArea.setFont(new Font("Monospaced", Font.PLAIN, 12));
            JScrollPane scrollPane = new JScrollPane(resultsArea);

            // Panneau pour les heatmaps
            JPanel heatmapsOuterPanel = new JPanel(new BorderLayout(5,5));
            JLabel heatmapTitleLabel = new JLabel("Heatmaps (Numérique, Analytique, Erreur)", SwingConstants.CENTER);
            heatmapsOuterPanel.add(heatmapTitleLabel, BorderLayout.NORTH);
            JPanel heatmapsGridPanel = new JPanel(new GridLayout(1,3,5,5)); // 1 ligne, 3 colonnes, avec espacement
            heatmapsOuterPanel.add(heatmapsGridPanel, BorderLayout.CENTER);

            // Split pane pour séparer texte et graphiques
            JSplitPane splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, scrollPane, heatmapsOuterPanel);
            splitPane.setResizeWeight(0.4); // Donne plus de place aux graphiques initialement

            mainFrame.add(controlPanel, BorderLayout.NORTH);
            mainFrame.add(splitPane, BorderLayout.CENTER);


            solveButton.addActionListener(e -> {
                UValueProvider2D selectedExactSolution = (UValueProvider2D) exactSolutionCombo.getSelectedItem();
                resultsArea.setText("");
                heatmapsGridPanel.removeAll();


                int[] N_values = {10, 20, 40, 80};
                List<Solution2D> solutions = new ArrayList<>();
                Solution2D prevSol = null;

                resultsArea.append("Résolution pour u_exact: " + selectedExactSolution.getName() + "\n");
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
                    heatmapsGridPanel.add(new HeatmapPanel(lastSol.numerical, "Numérique N=" + lastSol.nx));
                    heatmapsGridPanel.add(new HeatmapPanel(lastSol.analytical, "Analytique N=" + lastSol.nx));

                    double[][] errorGrid = new double[lastSol.nx+1][lastSol.ny+1];
                    for(int i=0; i<=lastSol.nx; i++) for(int j=0; j<=lastSol.ny; j++) errorGrid[i][j] = Math.abs(lastSol.numerical[i][j] - lastSol.analytical[i][j]);
                    heatmapsGridPanel.add(new HeatmapPanel(errorGrid, "Erreur Absolue N=" + lastSol.nx));
                }
                heatmapsGridPanel.revalidate();
                heatmapsGridPanel.repaint();
                mainFrame.pack();
            });
            
            mainFrame.setMinimumSize(new Dimension(600, 700)); // Taille minimale
            mainFrame.pack();
            mainFrame.setLocationRelativeTo(null);
            mainFrame.setVisible(true);
            
            // Lancer une résolution initiale pour exemple
            if (exactSolutionCombo.getItemCount() > 0) {
                 exactSolutionCombo.setSelectedIndex(0);
                 solveButton.doClick();
            }
        });
    }
}

class HeatmapPanel extends JPanel {
    private double[][] data;
    private String title;
    private double minVal = Double.MAX_VALUE, maxVal = Double.MIN_VALUE;

    public HeatmapPanel(double[][] dataGrid, String panelTitle) {
        this.data = dataGrid;
        this.title = panelTitle;
        setPreferredSize(new Dimension(250, 280)); // Taille réduite pour tenir dans GridLayout

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
        if (minVal == maxVal) maxVal = minVal + 1e-9; // S'assurer que maxVal > minVal
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        if (data == null || data.length == 0 || data[0].length == 0) {
            g.drawString("Pas de données", 10, 20);
            return;
        }

        int panelWidth = getWidth();
        int panelHeight = getHeight();
        // Laisser de la place pour titre et légende
        int usableWidth = panelWidth - 20; // 10px marge de chaque côté
        int usableHeight = panelHeight - 40; // 20px pour titre, 20px pour légende

        int rows = data.length;      // Correspond à ny+1 points (direction y de la grille)
        int cols = data[0].length;   // Correspond à nx+1 points (direction x de la grille)

        int cellWidth = Math.max(1, usableWidth / cols);
        int cellHeight = Math.max(1, usableHeight / rows);
        
        // Centre la heatmap
        int offsetX = (panelWidth - cols * cellWidth) / 2;
        int offsetY = 20 + (usableHeight - rows * cellHeight) / 2; // 20 pour titre

        g.setColor(Color.BLACK);
        g.drawString(title, panelWidth/2 - g.getFontMetrics().stringWidth(title)/2, 15);

        for (int i = 0; i < rows; i++) { // i itère sur les lignes de la matrice `data` (souvent direction y)
            for (int j = 0; j < cols; j++) { // j itère sur les colonnes (souvent direction x)
                double value = data[i][j]; // data[y_idx][x_idx] si on veut (i=y, j=x)
                float normalized = 0.5f;
                if (maxVal > minVal) { // Eviter division par zero
                     normalized = (float) ((value - minVal) / (maxVal - minVal));
                }
                normalized = Math.max(0f, Math.min(1f, normalized)); // Clamp entre 0 et 1

                Color color;
                if (normalized < 0.5f) {
                    color = new Color(normalized * 2, normalized * 2, 1f); // Bleu -> Cyan -> Blanc-ish
                } else {
                     color = new Color(1f, (1 - normalized) * 2, (1 - normalized) * 2); // Blanc-ish -> Jaune -> Rouge
                }
                g.setColor(color);
                // Affichage standard : data[i][j] où i est l'index de ligne (y), j est l'index de colonne (x)
                // L'origine du graphique (0,0) est en haut à gauche.
                // Pour que l'axe y pointe vers le haut, on inverse la coordonnée y du dessin.
                // (rows - 1 - i) pour inverser l'axe des y si data[0] est la ligne du bas.
                // Si data[0] est la ligne du haut (comme souvent en image), alors juste 'i'.
                // Ici, on suppose que data[i][j] correspond à u(x_j, y_i) où y_i est croissant vers le haut.
                // Donc, pour dessiner, la ligne i=0 (y=0) doit être en bas.
                g.fillRect(offsetX + j * cellWidth, offsetY + (rows - 1 - i) * cellHeight, cellWidth, cellHeight);
            }
        }
        g.setColor(Color.BLACK);
        g.drawString(String.format("Min: %.2e", minVal), offsetX, panelHeight - 5);
        g.drawString(String.format("Max: %.2e", maxVal), offsetX + cols * cellWidth - g.getFontMetrics().stringWidth(String.format("Max: %.2e", maxVal)), panelHeight - 5);
    }
}
