/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package FiniteDifference1;

/**
 *
 * @author Roddier
 */
import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;


public class ODEFiniteDifference {

    // Interface pour la solution exacte u(x)
    interface UValueProvider {
        double getValue(double x); // u(x)
        double getFirstDerivative(double x); // u'(x)
        double getSecondDerivative(double x); // u''(x)
        String getName(); // Nom de u(x), ex: "sin(πx)"
    }

    // Implémentations de UValueProvider
    static class USinPiX implements UValueProvider {
        public double getValue(double x) { return Math.sin(Math.PI * x); }
        public double getFirstDerivative(double x) { return Math.PI * Math.cos(Math.PI * x); }
        public double getSecondDerivative(double x) { return -Math.PI * Math.PI * Math.sin(Math.PI * x); }
        public String getName() { return "u(x) = sin(πx)"; }
    }

    static class UXCube implements UValueProvider {
        public double getValue(double x) { return x * x * x; }
        public double getFirstDerivative(double x) { return 3 * x * x; }
        public double getSecondDerivative(double x) { return 6 * x; }
        public String getName() { return "u(x) = x³"; }
    }
    
    // Types d'équations différentielles
    enum EquationType {
        TYPE1("-u'' + u = f (Source: -u_exact'' + u_exact)"),
        TYPE2("-u'' + u' = f (Source: -u_exact'' + u_exact')"),
        TYPE3("-u'' = f (Source: -u_exact'')");
        
        private final String description;
        EquationType(String description) { this.description = description; }
        public String getDescription() { return description; }
    }
    
    // Classe pour stocker les résultats
    static class Solution {
        double[] x;
        double[] numerical;
        double[] analytical;
        double error;
        int n;
        EquationType type;
        UValueProvider exactSolution;
        
        Solution(double[] x, double[] numerical, double[] analytical, double error, int n, EquationType type, UValueProvider exactSolution) {
            this.x = x;
            this.numerical = numerical;
            this.analytical = analytical;
            this.error = error;
            this.n = n;
            this.type = type;
            this.exactSolution = exactSolution;
        }
    }
    
    private static double[] solveTridiagonal(double[] a, double[] b, double[] c, double[] d) {
        int n_eq = b.length;
        if (n_eq == 0) return new double[0];
        double[] cp = new double[n_eq];
        double[] dp = new double[n_eq];
        double[] x_sol = new double[n_eq];
        
        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];
        
        for (int i = 1; i < n_eq; i++) {
            double denom = b[i] - a[i] * cp[i-1];
            cp[i] = c[i] / denom;
            dp[i] = (d[i] - a[i] * dp[i-1]) / denom;
        }
        
        x_sol[n_eq-1] = dp[n_eq-1];
        for (int i = n_eq-2; i >= 0; i--) {
            x_sol[i] = dp[i] - cp[i] * x_sol[i+1];
        }
        
        return x_sol;
    }
    
    private static double getAnalyticalSolutionValue(double x_coord, UValueProvider uExactSource, EquationType type, double u0, double u1) {
        if (type == EquationType.TYPE3) {
            double uExactValAt0 = uExactSource.getValue(0.0);
            double uExactValAt1 = uExactSource.getValue(1.0);
            double c0 = u0 - uExactValAt0;
            double c1 = u1 - uExactValAt1 - c0;
            return uExactSource.getValue(x_coord) + c1 * x_coord + c0;
        } else {
            System.err.println("Avertissement: la solution analytique pour " + type +
                               " suppose que u0/u1 correspondent aux valeurs de uExactSource aux bords.");
            return uExactSource.getValue(x_coord);
        }
    }
    
    public static Solution solve(int n_intervals, EquationType type, UValueProvider uExact, double u0_bc, double u1_bc) {
        double h_step = 1.0 / n_intervals;
        double[] x_nodes = new double[n_intervals+1];
        for (int i = 0; i <= n_intervals; i++) {
            x_nodes[i] = i * h_step;
        }
        
        if (n_intervals <= 1) {
            double[] numerical = new double[n_intervals+1];
            double[] analytical = new double[n_intervals+1];
            if (n_intervals==0) {
                 if (x_nodes.length > 0) {
                    numerical[0] = (u0_bc + u1_bc) / 2.0;
                    analytical[0] = getAnalyticalSolutionValue(x_nodes[0], uExact, type, u0_bc, u1_bc);
                 }
            } else {
                numerical[0] = u0_bc;
                numerical[1] = u1_bc;
            }
            for (int i = 0; i <= n_intervals; i++) {
                 if (i < x_nodes.length)
                    analytical[i] = getAnalyticalSolutionValue(x_nodes[i], uExact, type, u0_bc, u1_bc);
            }
            double max_abs_error_base_case = 0;
            for (int i = 0; i <= n_intervals; i++) {
                 if (i < numerical.length && i < analytical.length) {
                    double current_abs_error = Math.abs(numerical[i] - analytical[i]);
                    if (current_abs_error > max_abs_error_base_case) {
                        max_abs_error_base_case = current_abs_error;
                    }
                }
            }
            return new Solution(x_nodes, numerical, analytical, max_abs_error_base_case, n_intervals, type, uExact);
        }

        int n_unknowns = n_intervals - 1;
        double[] a_sub = new double[n_unknowns];
        double[] b_diag = new double[n_unknowns];
        double[] c_sur = new double[n_unknowns];
        double[] d_rhs = new double[n_unknowns];

        java.util.function.Function<Double, Double> f_provider;
        switch (type) {
            case TYPE1:
                f_provider = (x_val) -> -uExact.getSecondDerivative(x_val) + uExact.getValue(x_val);
                break;
            case TYPE2:
                f_provider = (x_val) -> -uExact.getSecondDerivative(x_val) + uExact.getFirstDerivative(x_val);
                break;
            case TYPE3:
            default:
                f_provider = (x_val) -> -uExact.getSecondDerivative(x_val);
                break;
        }

        switch (type) {
            case TYPE1:
                for (int i = 0; i < n_unknowns; i++) {
                    int actual_node_idx = i + 1;
                    a_sub[i] = -1.0/(h_step*h_step);
                    b_diag[i] = 2.0/(h_step*h_step) + 1.0;
                    c_sur[i] = -1.0/(h_step*h_step);
                    d_rhs[i] = f_provider.apply(x_nodes[actual_node_idx]);
                }
                d_rhs[0] -= a_sub[0] * u0_bc;
                if (n_unknowns > 0) a_sub[0] = 0;
                if (n_unknowns > 1) d_rhs[n_unknowns-1] -= c_sur[n_unknowns-1] * u1_bc;
                if (n_unknowns > 0) c_sur[n_unknowns-1] = 0;
                break;
            case TYPE2:
                for (int i = 0; i < n_unknowns; i++) {
                    int actual_node_idx = i + 1;
                    a_sub[i] = -1.0/(h_step*h_step) - 1.0/(2*h_step);
                    b_diag[i] = 2.0/(h_step*h_step);
                    c_sur[i] = -1.0/(h_step*h_step) + 1.0/(2*h_step);
                    d_rhs[i] = f_provider.apply(x_nodes[actual_node_idx]);
                }
                d_rhs[0] -= a_sub[0] * u0_bc;
                if (n_unknowns > 0) a_sub[0] = 0;
                if (n_unknowns > 1) d_rhs[n_unknowns-1] -= c_sur[n_unknowns-1] * u1_bc;
                if (n_unknowns > 0) c_sur[n_unknowns-1] = 0;
                break;
            case TYPE3:
                for (int i = 0; i < n_unknowns; i++) {
                    int actual_node_idx = i + 1;
                    a_sub[i] = -1.0;
                    b_diag[i] = 2.0;
                    c_sur[i] = -1.0;
                    d_rhs[i] = h_step*h_step * f_provider.apply(x_nodes[actual_node_idx]);
                }
                d_rhs[0] += u0_bc;
                if (n_unknowns > 0) a_sub[0] = 0;
                if (n_unknowns > 1) d_rhs[n_unknowns-1] += u1_bc;
                else if (n_unknowns == 1) d_rhs[0] += u1_bc;
                if (n_unknowns > 0) c_sur[n_unknowns-1] = 0;
                break;
        }
        
        double[] u_inner = solveTridiagonal(a_sub, b_diag, c_sur, d_rhs);
        
        double[] numerical = new double[n_intervals+1];
        double[] analytical = new double[n_intervals+1];
        
        numerical[0] = u0_bc;
        numerical[n_intervals] = u1_bc;
        for (int i = 0; i < u_inner.length; i++) {
            numerical[i+1] = u_inner[i];
        }
        
        for (int i = 0; i <= n_intervals; i++) {
            analytical[i] = getAnalyticalSolutionValue(x_nodes[i], uExact, type, u0_bc, u1_bc);
        }
        
        double max_abs_error = 0;
        for (int i = 0; i <= n_intervals; i++) {
            double current_abs_error = Math.abs(numerical[i] - analytical[i]);
            if (current_abs_error > max_abs_error) {
                max_abs_error = current_abs_error;
            }
        }
        
        return new Solution(x_nodes, numerical, analytical, max_abs_error, n_intervals, type, uExact);
    }
    
    public static double calculateConvergenceOrder(double error1, double error2, int n1_intervals, int n2_intervals) {
        if (n1_intervals == n2_intervals) {
            if (Math.abs(error1 - error2) < 1e-12) return Double.NaN;
            return (error1 > error2) ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
        }
        if (error1 <= 1e-14 && error2 <= 1e-14) return Double.NaN;
        if (error2 <= 1e-14) return Double.POSITIVE_INFINITY;
        if (error1 <= 1e-14) return Double.NEGATIVE_INFINITY;

        double ratio_err = error1 / error2;
        double ratio_n = (double)n2_intervals / n1_intervals;
        if (ratio_err <= 0 || ratio_n <=0) return Double.NaN;

        return Math.log(ratio_err) / Math.log(ratio_n);
    }
    
    static class GraphPanel extends JPanel {
        private List<Solution> solutions;
        private boolean showErrorGraph;
        
        public GraphPanel(List<Solution> solutions, boolean showErrorGraph) {
            this.solutions = solutions;
            this.showErrorGraph = showErrorGraph;
            setPreferredSize(new Dimension(800, 600));
            setBackground(Color.WHITE);
        }
        
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            
            int plotWidth = getWidth() - 100;
            int plotHeight = getHeight() - 100;
            int margin = 50;
            
            g2.setColor(Color.BLACK);
            g2.drawLine(margin, plotHeight + margin, plotWidth + margin, plotHeight + margin);
            g2.drawLine(margin, margin, margin, plotHeight + margin);
            
            if (showErrorGraph) {
                drawErrorGraph(g2, margin, plotWidth, plotHeight);
            } else {
                drawSolutionGraph(g2, margin, plotWidth, plotHeight);
            }
        }
        
        private void drawSolutionGraph(Graphics2D g2, int margin, int plotWidth, int plotHeight) {
            if (solutions.isEmpty()) return;
            
            Solution sol = solutions.get(0);
            
            double minY = Double.MAX_VALUE, maxY = Double.MIN_VALUE;
            for (int i = 0; i < sol.x.length; i++) {
                minY = Math.min(minY, Math.min(sol.numerical[i], sol.analytical[i]));
                maxY = Math.max(maxY, Math.max(sol.numerical[i], sol.analytical[i]));
            }
             if (Math.abs(maxY - minY) < 1e-9) {
                minY -= 0.5;
                maxY += 0.5;
            }
            double y_range = maxY - minY;
            
            g2.setColor(Color.BLUE);
            g2.setStroke(new BasicStroke(2));
            for (int i = 0; i < sol.x.length - 1; i++) {
                int x1 = margin + (int)(sol.x[i] * plotWidth);
                int y1 = margin + plotHeight - (int)((sol.analytical[i] - minY) / y_range * plotHeight);
                int x2 = margin + (int)(sol.x[i+1] * plotWidth);
                int y2 = margin + plotHeight - (int)((sol.analytical[i+1] - minY) / y_range * plotHeight);
                g2.drawLine(x1, y1, x2, y2);
            }
            
            g2.setColor(Color.RED);
            g2.setStroke(new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0, new float[]{5}, 0));
            for (int i = 0; i < sol.x.length - 1; i++) {
                int x1 = margin + (int)(sol.x[i] * plotWidth);
                int y1 = margin + plotHeight - (int)((sol.numerical[i] - minY) / y_range * plotHeight);
                int x2 = margin + (int)(sol.x[i+1] * plotWidth);
                int y2 = margin + plotHeight - (int)((sol.numerical[i+1] - minY) / y_range * plotHeight);
                g2.drawLine(x1, y1, x2, y2);
            }
            
            g2.setColor(Color.BLACK);
            g2.setFont(new Font("Arial", Font.BOLD, 16));
            g2.drawString(sol.type.getDescription() + ", Solution exacte u(x): " + sol.exactSolution.getName(), margin, 25); // "Exact u(x)" -> "Solution exacte u(x)"
            g2.drawString(String.format("N = %d, C.L.: u(0)=%.2f, u(1)=%.2f", sol.n, sol.numerical[0], sol.numerical[sol.n]), margin, 45); // BC -> C.L.
            
            g2.setFont(new Font("Arial", Font.PLAIN, 12));
            g2.setColor(Color.BLUE);
            g2.drawString("— Solution analytique", plotWidth + margin - 150, margin + 20); // "Analytical Solution" -> "Solution analytique"
            g2.setColor(Color.RED);
            g2.drawString("--- Solution numérique", plotWidth + margin - 150, margin + 40); // "Numerical Solution" -> "Solution numérique"
            g2.setColor(Color.BLACK);
            g2.drawString(String.format("Erreur L∞: %.2e", sol.error), plotWidth + margin - 150, margin + 60); // "L∞ Error" -> "Erreur L∞"
        }
        
        private void drawErrorGraph(Graphics2D g2, int margin, int plotWidth, int plotHeight) {
            if (solutions.size() < 2 && (solutions.isEmpty() || solutions.get(0).error <= 1e-18) ) {
                g2.drawString("Pas assez de données ou erreur trop faible pour le graphique d'erreur.", margin, plotHeight / 2); // Traduit
                return;
            }
             if (solutions.size() == 1 && solutions.get(0).error > 1e-18) {
                Solution sol = solutions.get(0);
                g2.setColor(Color.BLUE);
                int x_pt = margin + plotWidth / 2;
                int y_pt = margin + plotHeight / 2;
                g2.fillOval(x_pt - 3, y_pt - 3, 6, 6);
                g2.drawString(String.format("N=%d, Erreur=%.2e", sol.n, sol.error), x_pt + 5, y_pt + 5); // Error -> Erreur
                g2.setColor(Color.BLACK);
                g2.setFont(new Font("Arial", Font.BOLD, 16));
                g2.drawString("Erreur L∞ pour N=" + sol.n + " - " + sol.type.getDescription(), margin, 25); // L∞ Error -> Erreur L∞
                g2.drawString("Solution exacte u(x): " + sol.exactSolution.getName(), margin, 45); // Exact Solution -> Solution exacte
                return;
            }

            List<Solution> positiveErrorSolutions = solutions.stream()
                                                             .filter(s -> s != null && s.error > 1e-18)
                                                             .collect(Collectors.toList());
            if (positiveErrorSolutions.size() < 2) {
                 g2.drawString("Pas assez de points d'erreur positifs pour le graphique log-log.", margin, plotHeight / 2); // Traduit
                return;
            }

            double minError = positiveErrorSolutions.stream().mapToDouble(s -> s.error).min().orElse(1e-16);
            double maxError = positiveErrorSolutions.stream().mapToDouble(s -> s.error).max().orElse(1.0);
            
            if (Math.abs(maxError - minError) < 1e-9 * Math.abs(maxError)) {
                 minError = maxError / ( (maxError > 1e-15) ? 10 : 1e-1);
                 if (minError == 0 && maxError == 0) {minError = 1e-16; maxError = 1e-15;}
                 else if (minError == 0) minError = maxError * 1e-1;
                 if (minError > maxError) minError = maxError /10;
                 if (minError <=0) minError = 1e-16;
            }

            double logMinError = Math.log10(minError);
            double logMaxError = Math.log10(maxError);
            double logErrorRange = logMaxError - logMinError;
            if (Math.abs(logErrorRange) < 1e-9) logErrorRange = (logMinError != 0) ? Math.abs(logMinError * 0.1) : 1.0;
            
            int minN = positiveErrorSolutions.stream().mapToInt(s -> s.n).min().orElse(10);
            int maxN = positiveErrorSolutions.stream().mapToInt(s -> s.n).max().orElse(1000);
            double logMinN = Math.log10(minN);
            double logMaxN = Math.log10(maxN);
            double logNRange = logMaxN - logMinN;
            if (Math.abs(logNRange) < 1e-9) logNRange = (logMinN !=0) ? Math.abs(logMinN * 0.1) : 1.0;
            
            g2.setColor(Color.BLUE);
            g2.setStroke(new BasicStroke(2));
            
            Path2D.Double errorPath = new Path2D.Double();
            Solution firstPositiveErrorSol = positiveErrorSolutions.get(0);
            double firstX = margin + (int)((Math.log10(firstPositiveErrorSol.n) - logMinN) / logNRange * plotWidth);
            double firstY = margin + plotHeight - (int)((Math.log10(firstPositiveErrorSol.error) - logMinError) / logErrorRange * plotHeight);
            errorPath.moveTo(firstX,firstY);
            g2.fillOval((int)firstX - 3, (int)firstY - 3, 6, 6);

            for (int i = 1; i < positiveErrorSolutions.size(); i++) {
                Solution s = positiveErrorSolutions.get(i);
                int x_val = margin + (int)((Math.log10(s.n) - logMinN) / logNRange * plotWidth);
                int y_val = margin + plotHeight - (int)((Math.log10(s.error) - logMinError) / logErrorRange * plotHeight);
                errorPath.lineTo(x_val, y_val);
                g2.fillOval(x_val - 3, y_val - 3, 6, 6);
            }
            g2.draw(errorPath);
            
            g2.setColor(Color.BLACK);
            g2.setFont(new Font("Arial", Font.BOLD, 16));
            Solution firstSolForTitle = solutions.get(0);
            g2.drawString("Évolution de l'erreur L∞ - " + firstSolForTitle.type.getDescription(), margin, 25); // Traduit
            g2.drawString("Solution exacte u(x): " + firstSolForTitle.exactSolution.getName(), margin, 45); // Traduit
            
            g2.setFont(new Font("Arial", Font.PLAIN, 12));
            g2.drawString("Nombre d'intervalles N (échelle log)", plotWidth/2 + margin - 100, plotHeight + margin + 30); // Traduit
            
            AffineTransform orig = g2.getTransform();
            g2.rotate(-Math.PI/2);
            g2.drawString("Erreur L∞ (échelle log)", -(plotHeight/2 + margin + 50 ), margin - 30); // Traduit
            g2.setTransform(orig);
        }
    }
    
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            UValueProvider[] exactSolutions = {new USinPiX(), new UXCube()};
            EquationType[] types = {EquationType.TYPE3, EquationType.TYPE1, EquationType.TYPE2};
            
            JFrame mainFrame = new JFrame("Résolveur d'Équations Différentielles 1D (Différences Finies)"); // Traduit
            mainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            mainFrame.setLayout(new BorderLayout());
            
            JPanel controlPanel = new JPanel();
            JComboBox<UValueProvider> exactSolutionCombo = new JComboBox<>(exactSolutions);
            JComboBox<EquationType> typeCombo = new JComboBox<>(types);
            JButton solveButton = new JButton("Résoudre et Afficher Solution"); // Traduit
            JButton errorButton = new JButton("Afficher Courbe d'Erreur");    // Traduit
            
            exactSolutionCombo.setRenderer(new DefaultListCellRenderer() {
                @Override
                public Component getListCellRendererComponent(JList<?> list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
                    super.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);
                    if (value instanceof UValueProvider) {
                        setText(((UValueProvider) value).getName());
                    }
                    return this;
                }
            });
            typeCombo.setRenderer(new DefaultListCellRenderer() {
                @Override
                public Component getListCellRendererComponent(JList<?> list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
                    super.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);
                    if (value instanceof EquationType) {
                        setText(((EquationType) value).getDescription());
                    }
                    return this;
                }
            });

            controlPanel.add(new JLabel("Solution exacte u(x):")); // Traduit
            controlPanel.add(exactSolutionCombo);
            controlPanel.add(new JLabel("Type d'équation:"));    // Traduit
            controlPanel.add(typeCombo);
            controlPanel.add(solveButton);
            controlPanel.add(errorButton);
            
            JTabbedPane tabbedPane = new JTabbedPane();
            mainFrame.add(controlPanel, BorderLayout.NORTH);
            mainFrame.add(tabbedPane, BorderLayout.CENTER);
            
            final double u0_bc_val = 0.0;
            final double u1_bc_val = 1.0;

            solveButton.addActionListener(e -> {
                UValueProvider selectedExactSolution = (UValueProvider) exactSolutionCombo.getSelectedItem();
                EquationType selectedType = (EquationType) typeCombo.getSelectedItem();
                
                tabbedPane.removeAll();
                
                int[] meshSizesForDisplay = {10, 20, 40, 80};
                System.out.println("\n=== Solutions pour " + selectedType.getDescription() +
                                   " avec Solution Exacte u(x) = " + selectedExactSolution.getName() + " ==="); // Traduit
                System.out.println("Conditions aux Limites: u(0)=" + u0_bc_val + ", u(1)=" + u1_bc_val); // Traduit

                for (int n_val_loop : meshSizesForDisplay) {
                    Solution sol = solve(n_val_loop, selectedType, selectedExactSolution, u0_bc_val, u1_bc_val);
                    GraphPanel panel = new GraphPanel(Arrays.asList(sol), false);
                    tabbedPane.addTab("Solution (N = " + n_val_loop + ")", panel);
                    System.out.printf("N = %d: Erreur L∞ = %.6e\n", n_val_loop, sol.error); // Traduit
                }
            });

            errorButton.addActionListener(e -> {
                UValueProvider selectedExactSolution = (UValueProvider) exactSolutionCombo.getSelectedItem();
                EquationType selectedType = (EquationType) typeCombo.getSelectedItem();
                
                int[] meshSizesForErrorCurve = {10, 20, 40, 80, 160, 320};
                List<Solution> solutionsList = new ArrayList<>();
                
                System.out.println("\n=== Calcul de la Courbe d'Erreur pour " + selectedType.getDescription() +
                                   " avec Solution Exacte u(x) = " + selectedExactSolution.getName() + " ==="); // Traduit
                System.out.println("Conditions aux Limites: u(0)=" + u0_bc_val + ", u(1)=" + u1_bc_val); // Traduit

                Solution previousSolution = null;
                for (int n_val_loop : meshSizesForErrorCurve) {
                    Solution currentSolution = solve(n_val_loop, selectedType, selectedExactSolution, u0_bc_val, u1_bc_val);
                    solutionsList.add(currentSolution);

                    System.out.printf("N = %d: Erreur L∞ = %.6e\n", n_val_loop, currentSolution.error); // Traduit

                    if (previousSolution != null) {
                        double order = calculateConvergenceOrder(previousSolution.error, currentSolution.error,
                                                               previousSolution.n, currentSolution.n);
                        System.out.printf("Ordre de Convergence (entre N=%d et N=%d): %.2f\n", // Traduit
                                          previousSolution.n, currentSolution.n, order);
                    }
                    previousSolution = currentSolution;
                }
                
                System.out.println();
                
                GraphPanel errorCurvePanel = new GraphPanel(solutionsList, true);
                JFrame errorFrame = new JFrame("Évolution de l'Erreur L∞"); // Traduit
                errorFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
                errorFrame.add(errorCurvePanel);
                errorFrame.pack();
                errorFrame.setLocationRelativeTo(mainFrame);
                errorFrame.setVisible(true);
            });
            
            mainFrame.pack();
            mainFrame.setLocationRelativeTo(null);
            mainFrame.setVisible(true);
            
            exactSolutionCombo.setSelectedIndex(0);
            typeCombo.setSelectedIndex(0);
        });
    }
}