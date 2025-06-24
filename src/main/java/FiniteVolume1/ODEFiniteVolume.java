/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package FiniteVolume1; // Modifié

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


// Modifié ODEFiniteDifference en ODEFiniteVolume
public class ODEFiniteVolume {

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
        // u'(x) = πcos(πx), u''(x) = -π^2sin(πx)
        public double getSecondDerivative(double x) { return -Math.PI * Math.PI * Math.sin(Math.PI * x); }
        public String getName() { return "u(x) = sin(πx)"; }
    }

    static class UXCube implements UValueProvider {
        public double getValue(double x) { return x * x * x; }
        public double getFirstDerivative(double x) { return 3 * x * x; }
        // u'(x) = 3x^2, u''(x) = 6x
        public double getSecondDerivative(double x) { return 6 * x; }
        public String getName() { return "u(x) = x³"; }
    }
    
    // Types d'équations différentielles
    enum EquationType {
        TYPE1("-u'' + u = f"), // f = -u''_{exact} + u_{exact}
        TYPE2("-u'' + u' = f"), // f = -u''_{exact} + u'_{exact}
        TYPE3("-u'' = f");      // f = -u''_{exact}

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
        int n; // Nombre d'intervalles de maillage
        EquationType type;
        UValueProvider exactSolution; // Stocke la solution exacte u(x) utilisée
        
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
    
    // Résolution de système tridiagonal
    private static double[] solveTridiagonal(double[] a, double[] b, double[] c, double[] d) {
        int n_len = b.length;
        if (n_len == 0) return new double[0];
        double[] cp = new double[n_len];
        double[] dp = new double[n_len];
        double[] x_sol = new double[n_len];
        
        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];
        
        for (int i = 1; i < n_len; i++) {
            double denom = b[i] - a[i] * cp[i-1];
            cp[i] = c[i] / denom;
            dp[i] = (d[i] - a[i] * dp[i-1]) / denom;
        }
        
        x_sol[n_len-1] = dp[n_len-1];
        for (int i = n_len-2; i >= 0; i--) {
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
            System.err.println("Avertissement: getAnalyticalSolutionValue pour " + type + " suppose que u0/u1 correspondent à uExactSource aux bords.");
            return uExactSource.getValue(x_coord);
        }
    }
    
    public static Solution solve(int n_mesh, EquationType type, UValueProvider uExact, double u0_cl, double u1_cl) {
        double h = 1.0 / n_mesh;
        double[] x_coords = new double[n_mesh+1];
        for (int i = 0; i <= n_mesh; i++) {
            x_coords[i] = i * h;
        }
        
        if (n_mesh <= 1) {
            double[] numerical = new double[n_mesh+1];
            double[] analytical = new double[n_mesh+1];
            if (n_mesh==0) {
                 if (x_coords.length > 0) {
                    numerical[0] = (u0_cl+u1_cl)/2;
                    analytical[0] = getAnalyticalSolutionValue(x_coords[0], uExact, type, u0_cl, u1_cl);
                 }
            } else {
                numerical[0] = u0_cl;
                numerical[1] = u1_cl;
            }
            for (int i = 0; i <= n_mesh; i++) { // Boucle pour remplir analytical même pour n_mesh=0
                 if (i < analytical.length) // Protection supplémentaire
                    analytical[i] = getAnalyticalSolutionValue(x_coords[i], uExact, type, u0_cl, u1_cl);
            }
            double max_abs_error_base_case = 0;
            // Calcul d'erreur pour n_mesh=0 ou n_mesh=1
            for (int i = 0; i <= n_mesh; i++) {
                 if (i < numerical.length && i < analytical.length) { // Protection
                    double current_abs_error = Math.abs(numerical[i] - analytical[i]);
                    if (current_abs_error > max_abs_error_base_case) {
                        max_abs_error_base_case = current_abs_error;
                    }
                }
            }
            return new Solution(x_coords, numerical, analytical, max_abs_error_base_case, n_mesh, type, uExact);
        }

        double[] a_sub = new double[n_mesh-1];
        double[] b_diag = new double[n_mesh-1];
        double[] c_sur = new double[n_mesh-1];
        double[] d_rhs = new double[n_mesh-1];
        
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
                for (int i = 0; i < n_mesh-1; i++) {
                    int actual_idx = i + 1;
                    a_sub[i] = -1.0/(h*h);
                    b_diag[i] = 2.0/(h*h) + 1.0;
                    c_sur[i] = -1.0/(h*h);
                    d_rhs[i] = f_provider.apply(x_coords[actual_idx]);
                }
                d_rhs[0] -= a_sub[0] * u0_cl;
                if (n_mesh-1 > 0) a_sub[0] = 0;
                if (n_mesh-1 > 0) d_rhs[n_mesh-2] -= c_sur[n_mesh-2] * u1_cl;
                if (n_mesh-2 >=0 && n_mesh-2 < c_sur.length) c_sur[n_mesh-2] = 0;
                break;
            case TYPE2:
                for (int i = 0; i < n_mesh-1; i++) {
                    int actual_idx = i + 1;
                    a_sub[i] = -1.0/(h*h) - 1.0/(2*h);
                    b_diag[i] = 2.0/(h*h);
                    c_sur[i] = -1.0/(h*h) + 1.0/(2*h);
                    d_rhs[i] = f_provider.apply(x_coords[actual_idx]);
                }
                d_rhs[0] -= a_sub[0] * u0_cl;
                if (n_mesh-1 > 0) a_sub[0] = 0;
                if (n_mesh-1 > 0) d_rhs[n_mesh-2] -= c_sur[n_mesh-2] * u1_cl;
                if (n_mesh-2 >=0 && n_mesh-2 < c_sur.length) c_sur[n_mesh-2] = 0;
                break;
            case TYPE3:
                for (int i = 0; i < n_mesh-1; i++) {
                    int actual_idx = i + 1;
                    a_sub[i] = -1.0;
                    b_diag[i] = 2.0;
                    c_sur[i] = -1.0;
                    d_rhs[i] = h*h * f_provider.apply(x_coords[actual_idx]);
                }
                d_rhs[0] += u0_cl;
                if (n_mesh-1 > 0) a_sub[0] = 0;
                if (n_mesh-1 > 0) d_rhs[n_mesh-2] += u1_cl;
                if (n_mesh-2 >=0 && n_mesh-2 < c_sur.length) c_sur[n_mesh-2] = 0;
                break;
        }
        
        double[] u_inner = solveTridiagonal(a_sub, b_diag, c_sur, d_rhs);
        
        double[] numerical = new double[n_mesh+1];
        double[] analytical = new double[n_mesh+1];

        numerical[0] = u0_cl;
        numerical[n_mesh] = u1_cl;
        for (int i = 0; i < u_inner.length; i++) {
            numerical[i+1] = u_inner[i];
        }
        
        for (int i = 0; i <= n_mesh; i++) {
            analytical[i] = getAnalyticalSolutionValue(x_coords[i], uExact, type, u0_cl, u1_cl);
        }
        
        double max_abs_error = 0;
        for (int i = 0; i <= n_mesh; i++) {
            double current_abs_error = Math.abs(numerical[i] - analytical[i]);
            if (current_abs_error > max_abs_error) {
                max_abs_error = current_abs_error;
            }
        }

        return new Solution(x_coords, numerical, analytical, max_abs_error, n_mesh, type, uExact);
    }
    
    public static double calculateConvergenceOrder(double error1, double error2, int n1, int n2) {
        if (n1 == n2) { // Avoid division by zero in log((double)n2 / n1)
            if (Math.abs(error1 - error2) < 1e-12) return Double.NaN; // No change in N or error
            return (error1 > error2) ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY; // Error changed with same N
        }
        if (error1 <= 1e-12 && error2 <= 1e-12) return Double.NaN; // Both errors are effectively zero
        if (error2 <= 1e-12) return Double.POSITIVE_INFINITY; // Converged to zero
        if (error1 <= 1e-12) return Double.NEGATIVE_INFINITY; // Was zero, now has error (divergence from zero)

        double ratio_err = error1 / error2;
        double ratio_n = (double)n2 / n1;
        if (ratio_err <= 0 || ratio_n <=0) return Double.NaN; // Invalid log arguments

        return Math.log(ratio_err) / Math.log(ratio_n);
    }

    static class GraphPanel extends JPanel {
        private List<Solution> solutions;
        private boolean showError;
        
        public GraphPanel(List<Solution> solutions, boolean showError) {
            this.solutions = solutions;
            this.showError = showError;
            setPreferredSize(new Dimension(800, 600));
            setBackground(Color.WHITE);
        }
        
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            
            int width = getWidth() - 100;
            int height = getHeight() - 100;
            int margin = 50;
            
            g2.setColor(Color.BLACK);
            g2.drawLine(margin, height + margin, width + margin, height + margin);
            g2.drawLine(margin, margin, margin, height + margin);
            
            if (showError) {
                drawErrorGraph(g2, margin, width, height);
            } else {
                drawSolutionGraph(g2, margin, width, height);
            }
        }

        private void drawSolutionGraph(Graphics2D g2, int margin, int width, int height) {
            if (solutions.isEmpty() || solutions.get(0) == null || solutions.get(0).x == null || solutions.get(0).x.length == 0) return;
            
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

            double range = maxY - minY;
            
            g2.setColor(Color.BLUE);
            g2.setStroke(new BasicStroke(2));
            Path2D.Double analyticalPath = new Path2D.Double();
            analyticalPath.moveTo(margin + (sol.x[0] * width), margin + height - ((sol.analytical[0] - minY) / range * height));
            for (int i = 1; i < sol.x.length; i++) {
                 analyticalPath.lineTo(margin + (sol.x[i] * width), margin + height - ((sol.analytical[i] - minY) / range * height));
            }
            g2.draw(analyticalPath);
            
            g2.setColor(Color.RED);
            g2.setStroke(new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0, new float[]{5}, 0));
            Path2D.Double numericalPath = new Path2D.Double();
            numericalPath.moveTo(margin + (sol.x[0] * width), margin + height - ((sol.numerical[0] - minY) / range * height));
            for (int i = 1; i < sol.x.length; i++) {
                numericalPath.lineTo(margin + (sol.x[i] * width), margin + height - ((sol.numerical[i] - minY) / range * height));
            }
            g2.draw(numericalPath);
            
            g2.setColor(Color.BLACK);
            g2.setFont(new Font("Arial", Font.BOLD, 16));
            g2.drawString(sol.type.getDescription() + ", u_{exact}: " + sol.exactSolution.getName(), margin, 25);
            g2.drawString(String.format("N = %d, CL: u(0)=%.2f, u(1)=%.2f", sol.n, sol.numerical[0], sol.numerical[sol.n]), margin, 45);

            g2.setFont(new Font("Arial", Font.PLAIN, 12));
            g2.setColor(Color.BLUE);
            g2.drawString("— Solution analytique", width + margin - 150, margin + 20); // Ajusté position légende
            g2.setColor(Color.RED);
            g2.drawString("--- Solution numérique", width + margin - 150, margin + 40); // Ajusté position légende
            g2.setColor(Color.BLACK);
            g2.drawString(String.format("Erreur L∞: %.2e", sol.error), width + margin - 150, margin + 60); // Ajusté position légende
        }

        private void drawErrorGraph(Graphics2D g2, int margin, int width, int height) {
            if (solutions.isEmpty()) return;

            List<Solution> positiveErrorSolutions = solutions.stream()
                                                             .filter(s -> s != null && s.error > 1e-18) // Filtrer erreurs très petites ou nulles
                                                             .collect(Collectors.toList());
            if (positiveErrorSolutions.size() < 1) return; // Pas de points à dessiner

            if (positiveErrorSolutions.size() == 1) { // Afficher un seul point
                Solution sol = positiveErrorSolutions.get(0);
                g2.setColor(Color.BLUE);
                int x_pt = margin + width / 2;
                int y_pt = margin + height / 2;
                g2.fillOval(x_pt - 3, y_pt - 3, 6, 6);
                g2.drawString(String.format("N=%d, Err=%.2e", sol.n, sol.error), x_pt + 5, y_pt + 5);
                // Titre pour un seul point
                g2.setColor(Color.BLACK);
                g2.setFont(new Font("Arial", Font.BOLD, 16));
                g2.drawString("Erreur L∞ pour N=" + sol.n + " - " + sol.type.getDescription(), margin, 25);
                g2.drawString("u_{exact}: " + sol.exactSolution.getName(), margin, 45);
                return;
            }

            double minError = positiveErrorSolutions.stream().mapToDouble(s -> s.error).min().orElse(1e-16);
            double maxError = positiveErrorSolutions.stream().mapToDouble(s -> s.error).max().orElse(1.0);
            
            if (Math.abs(maxError - minError) < 1e-9 * Math.abs(maxError)) { // Si plage d'erreur relative très petite
                 minError = maxError / ( (maxError > 1e-15) ? 10 : 1e-1); // Créer une plage artificielle, éviter division par zéro
                 if (minError == 0 && maxError == 0) {minError = 1e-16; maxError = 1e-15;} // Cas où tout est à zero
                 else if (minError == 0) minError = maxError * 1e-1; // Si maxError non nul mais minError l'est devenu
                 if (minError > maxError) minError = maxError /10; // s'assurer min < max
                 if (minError <=0) minError = 1e-16; // S'assurer que minError est positif
            }


            double logMinError = Math.log10(minError);
            double logMaxError = Math.log10(maxError);
            double logRange = logMaxError - logMinError;
            if (Math.abs(logRange) < 1e-9) logRange = (logMinError != 0) ? Math.abs(logMinError * 0.1) : 1.0; // Eviter division par zero ou plage nulle

            
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
            double firstX = margin + (int)((Math.log10(firstPositiveErrorSol.n) - logMinN) / logNRange * width);
            double firstY = margin + height - (int)((Math.log10(firstPositiveErrorSol.error) - logMinError) / logRange * height);
            errorPath.moveTo(firstX,firstY);
            g2.fillOval((int)firstX - 3, (int)firstY - 3, 6, 6);

            for (int i = 1; i < positiveErrorSolutions.size(); i++) {
                Solution s = positiveErrorSolutions.get(i);
                int x_val = margin + (int)((Math.log10(s.n) - logMinN) / logNRange * width);
                int y_val = margin + height - (int)((Math.log10(s.error) - logMinError) / logRange * height);
                errorPath.lineTo(x_val, y_val);
                g2.fillOval(x_val - 3, y_val - 3, 6, 6);
            }
            g2.draw(errorPath);

            g2.setColor(Color.BLACK);
            g2.setFont(new Font("Arial", Font.BOLD, 16));
            Solution firstSolForTitle = solutions.get(0);
            g2.drawString("Évolution de l'erreur L∞ - " + firstSolForTitle.type.getDescription(), margin, 25);
            g2.drawString("Solution exacte u(x): " + firstSolForTitle.exactSolution.getName(), margin, 45);

            g2.setFont(new Font("Arial", Font.PLAIN, 12));
            g2.drawString("Nombre de mailles N (log)", width/2 + margin/2, height + margin + 30);

            AffineTransform orig = g2.getTransform();
            g2.rotate(-Math.PI/2);
            g2.drawString("Erreur L∞ (log)", -(height/2 + margin/2), margin - 10);
            g2.setTransform(orig);
        }
    }
    
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            UValueProvider[] exactSolutions = {new USinPiX(), new UXCube()};
            EquationType[] types = {EquationType.TYPE3, EquationType.TYPE1, EquationType.TYPE2}; // TYPE3 first for default
            
            JFrame mainFrame = new JFrame("1D Differential Equation Solver (Finite Volumes)"); // Translated
            mainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            mainFrame.setLayout(new BorderLayout());
            
            JPanel controlPanel = new JPanel();
            JComboBox<UValueProvider> exactSolutionCombo = new JComboBox<>(exactSolutions);
            JComboBox<EquationType> typeCombo = new JComboBox<>(types);
            JButton solveButton = new JButton("Solve & Display Solution"); // English
            JButton errorButton = new JButton("Display Error Curve");    // English
            
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

            controlPanel.add(new JLabel("Exact Solution u(x):")); // English
            controlPanel.add(exactSolutionCombo);
            controlPanel.add(new JLabel("Equation Type:"));    // English
            controlPanel.add(typeCombo);
            controlPanel.add(solveButton);
            controlPanel.add(errorButton);
            
            JTabbedPane tabbedPane = new JTabbedPane();
            mainFrame.add(controlPanel, BorderLayout.NORTH);
            mainFrame.add(tabbedPane, BorderLayout.CENTER);
            
            final double u0_cl_val = 0.0;
            final double u1_cl_val = 1.0;

            solveButton.addActionListener(e -> {
                UValueProvider selectedExactSolution = (UValueProvider) exactSolutionCombo.getSelectedItem();
                EquationType selectedType = (EquationType) typeCombo.getSelectedItem();

                tabbedPane.removeAll();

                int[] meshSizesForDisplay = {10, 20, 40, 80};
                System.out.println("\n=== Solutions (1D Finite Volumes) for " + selectedType.getDescription() +
                                   " with Exact u(x) = " + selectedExactSolution.getName() + " ==="); // Translated
                System.out.println("Boundary Conditions: u(0)=" + u0_cl_val + ", u(1)=" + u1_cl_val); // English

                for (int n_val_loop : meshSizesForDisplay) {
                    Solution sol = solve(n_val_loop, selectedType, selectedExactSolution, u0_cl_val, u1_cl_val);
                    GraphPanel panel = new GraphPanel(Arrays.asList(sol), false);
                    tabbedPane.addTab("Solution (N = " + n_val_loop + ")", panel); // Tab title can remain short
                    System.out.printf("N = %d: L∞ Error = %.6e\n", n_val_loop, sol.error); // English
                }
            });

            errorButton.addActionListener(e -> {
                UValueProvider selectedExactSolution = (UValueProvider) exactSolutionCombo.getSelectedItem();
                EquationType selectedType = (EquationType) typeCombo.getSelectedItem();

                tabbedPane.removeAll();

                int[] meshSizesForErrorCurve = {10, 20, 40, 80, 160, 320};
                List<Solution> solutionsList = new ArrayList<>(); // Renamed

                System.out.println("\n=== Calculating Error Curve (1D Finite Volumes) for " + selectedType.getDescription() +
                                   " with Exact u(x) = " + selectedExactSolution.getName() + " ==="); // Translated
                System.out.println("Boundary Conditions: u(0)=" + u0_cl_val + ", u(1)=" + u1_cl_val); // English

                Solution previousSolution = null; // Renamed
                for (int n_val_loop : meshSizesForErrorCurve) {
                    Solution currentSolution = solve(n_val_loop, selectedType, selectedExactSolution, u0_cl_val, u1_cl_val); // Renamed
                    solutionsList.add(currentSolution);
                    System.out.printf("N = %d: L∞ Error = %.6e\n", n_val_loop, currentSolution.error); // English
                    
                    if (previousSolution != null) {
                        double order = calculateConvergenceOrder(previousSolution.error, currentSolution.error,
                                                               previousSolution.n, currentSolution.n);
                        System.out.printf("Convergence Order (between N=%d and N=%d): %.2f\n", // English
                                          previousSolution.n, currentSolution.n, order);
                    }
                    previousSolution = currentSolution;
                }

                System.out.println();

                GraphPanel errorCurvePanel = new GraphPanel(solutionsList, true); // Renamed
                JFrame errorFrame = new JFrame("L∞ Error Evolution (1D Finite Volumes)"); // Translated
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
