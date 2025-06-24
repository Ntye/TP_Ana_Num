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

public class ODEFiniteDifference {

    // Interface for the exact solution u(x)
    interface UValueProvider {
        double getValue(double x); // u(x) value
        double getFirstDerivative(double x); // u'(x) value
        double getSecondDerivative(double x); // u''(x) value
        String getName(); // Name of u(x), e.g., "u(x) = sin(πx)"
    }

    // Implementations of UValueProvider
    static class USinPiX implements UValueProvider {
        public double getValue(double x) { return Math.sin(Math.PI * x); }
        public double getFirstDerivative(double x) { return Math.PI * Math.cos(Math.PI * x); }
        // u'(x) = πcos(πx), u''(x) = -π^2sin(πx)
        public double getSecondDerivative(double x) { return -Math.PI * Math.PI * Math.sin(Math.PI * x); }
        public String getName() { return "u(x) = sin(πx)"; } // Name used in UI
    }

    static class UXCube implements UValueProvider {
        public double getValue(double x) { return x * x * x; }
        public double getFirstDerivative(double x) { return 3 * x * x; }
        // u'(x) = 3x^2, u''(x) = 6x
        public double getSecondDerivative(double x) { return 6 * x; }
        public String getName() { return "u(x) = x³"; } // Name used in UI
    }
    
    // Types of differential equations
    enum EquationType {
        // Descriptions are mathematical forms, and also used in UI.
        TYPE1("-u'' + u = f (Source: -u_exact'' + u_exact)"),
        TYPE2("-u'' + u' = f (Source: -u_exact'' + u_exact')"),
        TYPE3("-u'' = f (Source: -u_exact'')");
        
        private final String description;
        EquationType(String description) { this.description = description; }
        public String getDescription() { return description; }
    }
    
    // Class to store results
    static class Solution {
        double[] x; // x-coordinates
        double[] numerical; // Numerical solution values
        double[] analytical; // Analytical solution values
        double error; // L-infinity error
        int n; // Number of mesh intervals
        EquationType type; // Type of equation solved
        UValueProvider exactSolution; // The exact solution u(x) used to derive f(x) for the RHS
        
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
    
    // Tridiagonal system solver (Thomas algorithm)
    private static double[] solveTridiagonal(double[] a, double[] b, double[] c, double[] d) {
        int n_eq = b.length; // n_eq: number of equations
        if (n_eq == 0) return new double[0]; // Handle empty system
        double[] cp = new double[n_eq]; // c' (c prime)
        double[] dp = new double[n_eq]; // d' (d prime)
        double[] x_sol = new double[n_eq]; // solution vector
        
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
    
    // Analytical Solutions
    // Calculates the analytical solution of the problem -u'' = f_src(x) with u(0)=u0, u(1)=u1,
    // where f_src(x) is derived from uExactSource (f_src(x) = -uExactSource.getSecondDerivative(x) for TYPE3).
    private static double getAnalyticalSolutionValue(double x_coord, UValueProvider uExactSource, EquationType type, double u0, double u1) {
        // u_particular(x) is an antiderivative of -f_src(x).
        // For TYPE3, -f_src(x) = uExactSource.getSecondDerivative(x).
        // So, a particular solution is uExactSource.getValue(x).
        // The general solution is u(x) = uExactSource.getValue(x) + C1*x + C0.

        // For TYPE3: -u'' = -uExactSource.getSecondDerivative(x)
        if (type == EquationType.TYPE3) {
            double uExactValAt0 = uExactSource.getValue(0.0);
            double uExactValAt1 = uExactSource.getValue(1.0);

            double c0 = u0 - uExactValAt0; // C0 = u0 - u_particular(0)
            double c1 = u1 - uExactValAt1 - c0; // C1 = u1 - u_particular(1) - C0

            return uExactSource.getValue(x_coord) + c1 * x_coord + c0;
        } else {
            // For TYPE1 and TYPE2, the derivation of f(x) is more complex if u0/u1 don't match uExactSource.
            // The current request focuses on -u''=f for specific uExact.
            // This part would need adjustment if other types were the primary focus with arbitrary BCs.
            System.err.println("Warning: Analytical solution for " + type +
                               " might not be accurate if u0/u1 don't match uExactSource boundary values.");
            return uExactSource.getValue(x_coord); // Fallback for non-TYPE3, assumes u0/u1 match uExactSource
        }
    }
    
    // Numerical solver method
    public static Solution solve(int n_intervals, EquationType type, UValueProvider uExact, double u0_bc, double u1_bc) { // n->n_intervals, u0/u1 -> u0_bc/u1_bc
        double h_step = 1.0 / n_intervals; // h -> h_step
        double[] x_nodes = new double[n_intervals+1]; // x -> x_nodes
        for (int i = 0; i <= n_intervals; i++) {
            x_nodes[i] = i * h_step;
        }
        
        // System matrices for Ax = b for the N-1 interior points u_1, ..., u_{N-1}
        // Array sizes will be n_intervals-1 if n_intervals > 1.
        if (n_intervals <= 1) { // Handles cases with no interior points or a single segment
            double[] numerical = new double[n_intervals+1];
            double[] analytical = new double[n_intervals+1];
            if (n_intervals==0) {
                 // Edge case: if n_intervals is 0, x_nodes has 1 point.
                 // Solution is ill-defined or could be an average. For now, assume n_intervals >= 1 for meaningful solve.
                 if (x_nodes.length > 0) { // Should be true if n_intervals=0 due to n_intervals+1 size
                    numerical[0] = (u0_bc + u1_bc) / 2.0; // Or some other convention for n=0
                    analytical[0] = getAnalyticalSolutionValue(x_nodes[0], uExact, type, u0_bc, u1_bc);
                 }
            } else { // n_intervals=1, two points x0, x1
                numerical[0] = u0_bc;
                numerical[1] = u1_bc;
            }
            // Calculate analytical solution at nodes
            for (int i = 0; i <= n_intervals; i++) {
                 if (i < x_nodes.length) // Ensure index is within bounds
                    analytical[i] = getAnalyticalSolutionValue(x_nodes[i], uExact, type, u0_bc, u1_bc);
            }
            // L-infinity error calculation for base cases
            double max_abs_error_base_case = 0;
            for (int i = 0; i <= n_intervals; i++) {
                 if (i < numerical.length && i < analytical.length) { // Ensure index is within bounds
                    double current_abs_error = Math.abs(numerical[i] - analytical[i]);
                    if (current_abs_error > max_abs_error_base_case) {
                        max_abs_error_base_case = current_abs_error;
                    }
                }
            }
            return new Solution(x_nodes, numerical, analytical, max_abs_error_base_case, n_intervals, type, uExact);
        }

        // Number of unknowns for interior points
        int n_unknowns = n_intervals - 1;
        double[] a_sub = new double[n_unknowns]; // sub-diagonal
        double[] b_diag = new double[n_unknowns]; // diagonal
        double[] c_sur = new double[n_unknowns]; // super-diagonal
        double[] d_rhs = new double[n_unknowns]; // right-hand side vector

        // Define f(x) for the right-hand side based on uExact and equation type
        java.util.function.Function<Double, Double> f_provider;
        switch (type) {
            case TYPE1: // f = -u_exact'' + u_exact
                f_provider = (x_val) -> -uExact.getSecondDerivative(x_val) + uExact.getValue(x_val);
                break;
            case TYPE2: // f = -u_exact'' + u_exact'
                f_provider = (x_val) -> -uExact.getSecondDerivative(x_val) + uExact.getFirstDerivative(x_val);
                break;
            case TYPE3: // f = -u_exact''
            default:    // Default to ensure f_provider is initialized
                f_provider = (x_val) -> -uExact.getSecondDerivative(x_val);
                break;
        }

        // Construct the system Ax = d_rhs
        switch (type) {
            case TYPE1: // -u'' + u = f  => (-1/h^2)u_i-1 + (2/h^2 + 1)u_i + (-1/h^2)u_i+1 = f_i
                for (int i = 0; i < n_unknowns; i++) {
                    int actual_node_idx = i + 1; // Index in the full x_nodes array
                    a_sub[i] = -1.0/(h_step*h_step);
                    b_diag[i] = 2.0/(h_step*h_step) + 1.0;
                    c_sur[i] = -1.0/(h_step*h_step);
                    d_rhs[i] = f_provider.apply(x_nodes[actual_node_idx]);
                }
                // Adjust for boundary conditions u0_bc and u1_bc
                d_rhs[0] -= a_sub[0] * u0_bc;
                if (n_unknowns > 0) a_sub[0] = 0; // Not used by Thomas for first row, but good practice

                if (n_unknowns > 1) d_rhs[n_unknowns-1] -= c_sur[n_unknowns-1] * u1_bc; // Check index for c_sur
                if (n_unknowns > 0) c_sur[n_unknowns-1] = 0; // Not used by Thomas for last row
                break;
                
            case TYPE2: // -u'' + u' = f => (-1/h^2-1/(2h))u_i-1 + (2/h^2)u_i + (-1/h^2+1/(2h))u_i+1 = f_i
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
                
            case TYPE3: // -u'' = f  =>  -u_i-1 + 2u_i - u_i+1 = h^2 * f_i
                for (int i = 0; i < n_unknowns; i++) {
                    int actual_node_idx = i + 1;
                    a_sub[i] = -1.0;
                    b_diag[i] = 2.0;
                    c_sur[i] = -1.0;
                    d_rhs[i] = h_step*h_step * f_provider.apply(x_nodes[actual_node_idx]);
                }
                // Adjust for boundary conditions
                d_rhs[0] += u0_bc;
                if (n_unknowns > 0) a_sub[0] = 0;

                if (n_unknowns > 1) d_rhs[n_unknowns-1] += u1_bc; // If only one unknown, u1 affects d_rhs[0]
                else if (n_unknowns == 1) d_rhs[0] += u1_bc; // Special case for N=2 (1 unknown)

                if (n_unknowns > 0) c_sur[n_unknowns-1] = 0;
                break;
        }
        
        // Solve the tridiagonal system for interior points
        double[] u_inner = solveTridiagonal(a_sub, b_diag, c_sur, d_rhs);
        
        // Construct the full numerical solution
        double[] numerical = new double[n_intervals+1];
        double[] analytical = new double[n_intervals+1];
        
        numerical[0] = u0_bc;
        numerical[n_intervals] = u1_bc;
        for (int i = 0; i < u_inner.length; i++) {
            numerical[i+1] = u_inner[i];
        }
        
        // Populate analytical solution array
        for (int i = 0; i <= n_intervals; i++) {
            analytical[i] = getAnalyticalSolutionValue(x_nodes[i], uExact, type, u0_bc, u1_bc);
        }
        
        // Calculate L-infinity error
        double max_abs_error = 0;
        for (int i = 0; i <= n_intervals; i++) {
            double current_abs_error = Math.abs(numerical[i] - analytical[i]);
            if (current_abs_error > max_abs_error) {
                max_abs_error = current_abs_error;
            }
        }
        
        return new Solution(x_nodes, numerical, analytical, max_abs_error, n_intervals, type, uExact);
    }
    
    // Calculation of numerical convergence order
    public static double calculateConvergenceOrder(double error1, double error2, int n1_intervals, int n2_intervals) { // n1/n2 -> n1_intervals/n2_intervals
        if (n1_intervals == n2_intervals) {
            if (Math.abs(error1 - error2) < 1e-12) return Double.NaN; // No change in N or error
            return (error1 > error2) ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY; // Error changed with same N
        }
        if (error1 <= 1e-14 && error2 <= 1e-14) return Double.NaN; // Both errors effectively zero
        if (error2 <= 1e-14) return Double.POSITIVE_INFINITY; // Converged to zero error
        if (error1 <= 1e-14) return Double.NEGATIVE_INFINITY; // Was zero, now has error (divergence from zero)

        double ratio_err = error1 / error2;
        double ratio_n = (double)n2_intervals / n1_intervals;
        if (ratio_err <= 0 || ratio_n <=0) return Double.NaN; // Invalid log arguments

        return Math.log(ratio_err) / Math.log(ratio_n);
    }
    
    // Class for graphical visualization
    // Class for graphical visualization
    static class GraphPanel extends JPanel {
        private List<Solution> solutions;
        private boolean showErrorGraph; // Renamed from showError to showErrorGraph for clarity
        
        public GraphPanel(List<Solution> solutions, boolean showErrorGraph) { // Parameter name updated
            this.solutions = solutions;
            this.showErrorGraph = showErrorGraph; // Updated assignment
            setPreferredSize(new Dimension(800, 600));
            setBackground(Color.WHITE);
        }
        
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            
            int plotWidth = getWidth() - 100; // Renamed for clarity
            int plotHeight = getHeight() - 100; // Renamed for clarity
            int margin = 50;
            
            // Draw Axes
            g2.setColor(Color.BLACK);
            g2.drawLine(margin, plotHeight + margin, plotWidth + margin, plotHeight + margin); // x-axis
            g2.drawLine(margin, margin, margin, plotHeight + margin); // y-axis
            
            if (showErrorGraph) { // Condition updated
                drawErrorGraph(g2, margin, plotWidth, plotHeight);
            } else {
                drawSolutionGraph(g2, margin, plotWidth, plotHeight);
            }
        }
        
        private void drawSolutionGraph(Graphics2D g2, int margin, int plotWidth, int plotHeight) { // Parameters updated
            if (solutions.isEmpty()) return;
            
            Solution sol = solutions.get(0);
            
            // Trouver les limites
            double minY = Double.MAX_VALUE, maxY = Double.MIN_VALUE;
            for (int i = 0; i < sol.x.length; i++) {
                minY = Math.min(minY, Math.min(sol.numerical[i], sol.analytical[i]));
                maxY = Math.max(maxY, Math.max(sol.numerical[i], sol.analytical[i]));
            }
            
            double range = maxY - minY;
            if (range == 0) range = 1;
            
            // Dessiner la solution analytique
            g2.setColor(Color.BLUE);
            g2.setStroke(new BasicStroke(2));
            for (int i = 0; i < sol.x.length - 1; i++) {
                int x1 = margin + (int)(sol.x[i] * width);
                int y1 = margin + height - (int)((sol.analytical[i] - minY) / range * height);
                int x2 = margin + (int)(sol.x[i+1] * width);
                int y2 = margin + height - (int)((sol.analytical[i+1] - minY) / range * height);
                g2.drawLine(x1, y1, x2, y2);
            }
            
            // Dessiner la solution numérique
            g2.setColor(Color.RED);
            g2.setStroke(new BasicStroke(1, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 0, new float[]{5}, 0));
            for (int i = 0; i < sol.x.length - 1; i++) {
                int x1 = margin + (int)(sol.x[i] * width);
                int y1 = margin + height - (int)((sol.numerical[i] - minY) / range * height);
                int x2 = margin + (int)(sol.x[i+1] * width);
                int y2 = margin + height - (int)((sol.numerical[i+1] - minY) / range * height);
                g2.drawLine(x1, y1, x2, y2);
            }
            
            // Legend and Info
            g2.setColor(Color.BLACK);
            g2.setFont(new Font("Arial", Font.BOLD, 16));
            g2.drawString(sol.type.getDescription() + ", Exact u(x): " + sol.exactSolution.getName(), margin, 25);
            g2.drawString(String.format("N = %d, BC: u(0)=%.2f, u(1)=%.2f", sol.n, sol.numerical[0], sol.numerical[sol.n]), margin, 45); // CL -> BC, .1f -> .2f
            
            g2.setFont(new Font("Arial", Font.PLAIN, 12));
            g2.setColor(Color.BLUE);
            g2.drawString("— Analytical Solution", width - 150, margin + 20);
            g2.setColor(Color.RED);
            g2.drawString("--- Numerical Solution", width - 150, margin + 40);
            g2.setColor(Color.BLACK);
            g2.drawString(String.format("L∞ Error: %.2e", sol.error), width - 150, margin + 60);
        }
        
        private void drawErrorGraph(Graphics2D g2, int margin, int width, int height) { // Renamed plotWidth/Height to width/height for consistency
            if (solutions.size() < 2 && (solutions.isEmpty() || solutions.get(0).error <= 1e-18) ) {
                g2.drawString("Not enough data points or error too small for error graph.", margin, height / 2);
                return;
            }
             if (solutions.size() == 1 && solutions.get(0).error > 1e-18) {
                Solution sol = solutions.get(0);
                g2.setColor(Color.BLUE);
                int x_pt = margin + width / 2;
                int y_pt = margin + height / 2;
                g2.fillOval(x_pt - 3, y_pt - 3, 6, 6);
                g2.drawString(String.format("N=%d, Error=%.2e", sol.n, sol.error), x_pt + 5, y_pt + 5);
                g2.setColor(Color.BLACK);
                g2.setFont(new Font("Arial", Font.BOLD, 16));
                g2.drawString("L∞ Error for N=" + sol.n + " - " + sol.type.getDescription(), margin, 25);
                g2.drawString("Exact Solution u(x): " + sol.exactSolution.getName(), margin, 45);
                return;
            }

            List<Solution> positiveErrorSolutions = solutions.stream()
                                                             .filter(s -> s != null && s.error > 1e-18)
                                                             .collect(Collectors.toList());
            if (positiveErrorSolutions.size() < 2) {
                 g2.drawString("Not enough positive error data points for log-log graph.", margin, height / 2);
                return;
            }
            
            // Filter for positive errors for log scale
            double minError = solutions.stream()
                                     .mapToDouble(s -> s.error)
                                     .filter(e -> e > 0) // Important pour l'échelle log
                                     .min().orElse(1e-16); // Petite valeur si toutes les erreurs sont nulles
            double maxError = solutions.stream()
                                     .mapToDouble(s -> s.error)
                                     .max().orElse(1.0);
            if (minError <= 0) minError = 1e-16; // Assurer une valeur positive pour log
            if (maxError <= 0) maxError = 1.0; // Assurer une valeur positive si minError était aussi 0
            if (maxError < minError) maxError = minError * 10; // Assurer maxError > minError
            
            // Échelle logarithmique
            double logMinError = Math.log10(minError);
            double logMaxError = Math.log10(maxError);
            double logRange = logMaxError - logMinError;
            if (logRange == 0) logRange = 1;
            
            int minN = solutions.stream().mapToInt(s -> s.n).min().orElse(10);
            int maxN = solutions.stream().mapToInt(s -> s.n).max().orElse(1000);
            double logMinN = Math.log10(minN);
            double logMaxN = Math.log10(maxN);
            double logNRange = logMaxN - logMinN;
            
            // Dessiner les points et les relier
            g2.setColor(Color.BLUE);
            g2.setStroke(new BasicStroke(2));
            
            for (int i = 0; i < solutions.size() - 1; i++) {
                Solution s1 = solutions.get(i);
                Solution s2 = solutions.get(i + 1);
                
                int x1 = margin + (int)((Math.log10(s1.n) - logMinN) / logNRange * width);
                int y1 = margin + height - (int)((Math.log10(s1.error) - logMinError) / logRange * height);
                int x2 = margin + (int)((Math.log10(s2.n) - logMinN) / logNRange * width);
                int y2 = margin + height - (int)((Math.log10(s2.error) - logMinError) / logRange * height);
                
                g2.drawLine(x1, y1, x2, y2);
                g2.fillOval(x1 - 3, y1 - 3, 6, 6);
            }
            
            // Dernier point
            Solution lastSol = solutions.get(solutions.size() - 1);
            int lastX = margin + (int)((Math.log10(lastSol.n) - logMinN) / logNRange * width);
            int lastY = margin + height - (int)((Math.log10(lastSol.error) - logMinError) / logRange * height);
            g2.fillOval(lastX - 3, lastY - 3, 6, 6);
            
            // Title
            g2.setColor(Color.BLACK);
            g2.setFont(new Font("Arial", Font.BOLD, 16));
            Solution firstSolForTitle = solutions.get(0); // Use the first solution for consistent title info
            g2.drawString("L∞ Error Evolution - " + firstSolForTitle.type.getDescription(), margin, 25);
            g2.drawString("Exact Solution u(x): " + firstSolForTitle.exactSolution.getName(), margin, 45);
            
            // Axis Labels
            g2.setFont(new Font("Arial", Font.PLAIN, 12));
            g2.drawString("Number of Intervals N (log scale)", width/2 + margin - 100, height + margin + 30); // Adjusted position
            
            // Rotate for Y-axis label
            AffineTransform orig = g2.getTransform();
            g2.rotate(-Math.PI/2);
            g2.drawString("L∞ Error (log scale)", -(height/2 + margin + 50 ), margin - 30); // Adjusted position
            g2.setTransform(orig);
        }
    }
    
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            UValueProvider[] exactSolutions = {new USinPiX(), new UXCube()};
            EquationType[] types = {EquationType.TYPE3, EquationType.TYPE1, EquationType.TYPE2}; // TYPE3 first for default
            
            JFrame mainFrame = new JFrame("1D Differential Equation Solver (Finite Differences)");
            mainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            mainFrame.setLayout(new BorderLayout());
            
            JPanel controlPanel = new JPanel();
            JComboBox<UValueProvider> exactSolutionCombo = new JComboBox<>(exactSolutions);
            JComboBox<EquationType> typeCombo = new JComboBox<>(types);
            JButton solveButton = new JButton("Solve & Display Solution");
            JButton errorButton = new JButton("Display Error Curve");
            
            exactSolutionCombo.setRenderer(new DefaultListCellRenderer() {
                @Override
                public Component getListCellRendererComponent(JList<?> list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
                    super.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);
                    if (value instanceof UValueProvider) {
                        setText(((UValueProvider) value).getName()); // Assumes getName() is suitable for display
                    }
                    return this;
                }
            });
            typeCombo.setRenderer(new DefaultListCellRenderer() {
                @Override
                public Component getListCellRendererComponent(JList<?> list, Object value, int index, boolean isSelected, boolean cellHasFocus) {
                    super.getListCellRendererComponent(list, value, index, isSelected, cellHasFocus);
                    if (value instanceof EquationType) {
                        setText(((EquationType) value).getDescription()); // Assumes getDescription() is suitable
                    }
                    return this;
                }
            });

            controlPanel.add(new JLabel("Exact Solution u(x):"));
            controlPanel.add(exactSolutionCombo);
            controlPanel.add(new JLabel("Equation Type:"));
            controlPanel.add(typeCombo);
            controlPanel.add(solveButton);
            controlPanel.add(errorButton);
            
            JTabbedPane tabbedPane = new JTabbedPane();
            mainFrame.add(controlPanel, BorderLayout.NORTH);
            mainFrame.add(tabbedPane, BorderLayout.CENTER);
            
            final double u0_bc_val = 0.0; // Renamed for clarity
            final double u1_bc_val = 1.0;

            solveButton.addActionListener(e -> {
                UValueProvider selectedExactSolution = (UValueProvider) exactSolutionCombo.getSelectedItem();
                EquationType selectedType = (EquationType) typeCombo.getSelectedItem();
                
                tabbedPane.removeAll();
                
                int[] meshSizesForDisplay = {10, 20, 40, 80};
                System.out.println("\n=== Solutions for " + selectedType.getDescription() +
                                   " with Exact u(x) = " + selectedExactSolution.getName() + " ==="); // u_exact -> Exact u(x)
                System.out.println("Boundary Conditions: u(0)=" + u0_bc_val + ", u(1)=" + u1_bc_val);

                for (int n_val_loop : meshSizesForDisplay) {
                    Solution sol = solve(n_val_loop, selectedType, selectedExactSolution, u0_bc_val, u1_bc_val);
                    GraphPanel panel = new GraphPanel(Arrays.asList(sol), false); // showErrorGraph = false
                    tabbedPane.addTab("Solution (N = " + n_val_loop + ")", panel);
                    System.out.printf("N = %d: L∞ Error = %.6e\n", n_val_loop, sol.error);
                }
            });

            errorButton.addActionListener(e -> {
                UValueProvider selectedExactSolution = (UValueProvider) exactSolutionCombo.getSelectedItem();
                EquationType selectedType = (EquationType) typeCombo.getSelectedItem();
                
                // Results for error curve calculation
                int[] meshSizesForErrorCurve = {10, 20, 40, 80, 160, 320};
                List<Solution> solutionsList = new ArrayList<>();

                System.out.println("\n=== Calculating Error Curve for " + selectedType.getDescription() +
                                   " with Exact u(x) = " + selectedExactSolution.getName() + " ==="); // u_exact -> Exact u(x)
                System.out.println("Boundary Conditions: u(0)=" + u0_bc_val + ", u(1)=" + u1_bc_val);
                
                Solution previousSolution = null;
                for (int n_val_loop : meshSizesForErrorCurve) {
                    Solution currentSolution = solve(n_val_loop, selectedType, selectedExactSolution, u0_bc_val, u1_bc_val);
                    solutionsList.add(currentSolution);

                    System.out.printf("N = %d: L∞ Error = %.6e\n", n_val_loop, currentSolution.error);

                    if (previousSolution != null) {
                        double order = calculateConvergenceOrder(previousSolution.error, currentSolution.error,
                                                               previousSolution.n, currentSolution.n);
                        System.out.printf("Convergence Order (between N=%d and N=%d): %.2f\n",
                                          previousSolution.n, currentSolution.n, order);
                    }
                    previousSolution = currentSolution;
                }
                
                System.out.println(); // Newline after table
                
                GraphPanel errorCurvePanel = new GraphPanel(solutionsList, true);
                JFrame errorFrame = new JFrame("L∞ Error Evolution");
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