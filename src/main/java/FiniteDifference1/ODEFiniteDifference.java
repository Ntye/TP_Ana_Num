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
        int n;
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
        int n = b.length;
        double[] cp = new double[n];
        double[] dp = new double[n];
        double[] x = new double[n];
        
        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];
        
        for (int i = 1; i < n; i++) {
            double denom = b[i] - a[i] * cp[i-1];
            cp[i] = c[i] / denom;
            dp[i] = (d[i] - a[i] * dp[i-1]) / denom;
        }
        
        x[n-1] = dp[n-1];
        for (int i = n-2; i >= 0; i--) {
            x[i] = dp[i] - cp[i] * x[i+1];
        }
        
        return x;
    }
    
    // Solutions analytiques
    // Calcule la solution analytique du problème -u'' = f_src(x) avec u(0)=u0, u(1)=u1,
    // où f_src(x) est dérivée de uExact (f_src(x) = -uExact.getSecondDerivative(x) pour TYPE3).
    private static double getAnalyticalSolutionValue(double x_coord, UValueProvider uExactSource, EquationType type, double u0, double u1) {
        // u_particular(x) est une primitive de -f_src(x).
        // Pour TYPE3, -f_src(x) = uExactSource.getSecondDerivative(x).
        // Donc, une solution particulière est uExactSource.getValue(x).
        // La solution générale est u(x) = uExactSource.getValue(x) + C1*x + C0.

        // Pour TYPE3: -u'' = -uExactSource.getSecondDerivative(x)
        if (type == EquationType.TYPE3) {
            double uExactValAt0 = uExactSource.getValue(0.0);
            double uExactValAt1 = uExactSource.getValue(1.0);

            // C0 = u0 - uExactValAt0
            double c0 = u0 - uExactValAt0;
            // C1 = u1 - uExactValAt1 - C0
            double c1 = u1 - uExactValAt1 - c0;

            return uExactSource.getValue(x_coord) + c1 * x_coord + c0;
        } else {
            // Pour TYPE1 et TYPE2, la dérivation de f(x) est plus complexe.
            // La demande se concentre sur -u''=f. Provisoire:
            // Si on voulait gérer cela génériquement, il faudrait que UValueProvider fournisse la primitive double de f.
            // Pour l'instant, si ce n'est pas TYPE3, on retourne la valeur brute de uExactSource comme avant,
            // ce qui serait correct si u0 et u1 correspondaient aux valeurs de uExactSource aux bords.
            // Cela nécessitera une révision si TYPE1/TYPE2 sont utilisés avec des CLs arbitraires.
            // Pour la tâche actuelle, ceci est suffisant.
            System.err.println("Avertissement: getAnalyticalSolutionValue pour " + type + " suppose que u0/u1 correspondent à uExactSource aux bords.");
            return uExactSource.getValue(x_coord);
        }
    }
    
    // Résolution numérique
    public static Solution solve(int n, EquationType type, UValueProvider uExact, double u0, double u1) {
        double h = 1.0 / n;
        double[] x = new double[n+1];
        for (int i = 0; i <= n; i++) {
            x[i] = i * h;
        }
        
        // Matrices pour le système Ax = b pour les N-1 points intérieurs u_1, ..., u_{N-1}
        // La taille des tableaux a, b, c, d sera n-1 si n > 1. Si n=1, pas de points intérieurs.
        if (n <= 1) { // Cas où il n'y a pas de points intérieurs ou un seul segment
            double[] numerical = new double[n+1];
            double[] analytical = new double[n+1];
            if (n==0) {
                 // Cas très simple, la solution est juste les CLs si n=0 n'a pas de sens ici.
                 // Supposons n>=1 pour la discrétisation.
            } else { // n=1, deux points x0, x1
                numerical[0] = u0; // u(x_coords[0])
                numerical[1] = u1; // u(x_coords[1])
            }
            // La solution analytique est simplement uExact aux points
            for (int i = 0; i <= n; i++) {
                 analytical[i] = getAnalyticalSolutionValue(x[i], uExact, type, u0, u1); // x_coords -> x
            }
            // Calcul de l'erreur L-infini
            double max_abs_error_base_case = 0;
            if (n > 0) {
                for (int i = 0; i <= n; i++) {
                    double current_abs_error = Math.abs(numerical[i] - analytical[i]);
                    if (current_abs_error > max_abs_error_base_case) {
                        max_abs_error_base_case = current_abs_error;
                    }
                }
            } else { // n=0, un seul point, x[0]. Erreur est |u0 - uExact(x[0])|
                 if (n == 0 && x.length > 0) { // Protection // x_coords -> x
                    max_abs_error_base_case = Math.abs(u0 - analytical[0]);
                 }
            }
            return new Solution(x, numerical, analytical, max_abs_error_base_case, n, type, uExact); // x_coords -> x
        }

        double[] a_sub = new double[n-1]; // sous-diagonale
        double[] b_diag = new double[n-1]; // diagonale
        double[] c_sur = new double[n-1]; // sur-diagonale (renommé pour clarté)
        double[] d_rhs = new double[n-1]; // second membre
        
        // Définition de f(x) pour le membre de droite, basée sur uExact et le type d'équation
        java.util.function.Function<Double, Double> f_provider;
        switch (type) {
            case TYPE1: // f = -u''_{exact} + u_{exact}
                f_provider = (x_val) -> -uExact.getSecondDerivative(x_val) + uExact.getValue(x_val);
                break;
            case TYPE2: // f = -u''_{exact} + u'_{exact}
                f_provider = (x_val) -> -uExact.getSecondDerivative(x_val) + uExact.getFirstDerivative(x_val);
                break;
            case TYPE3: // f = -u''_{exact}
            default:    // Cas par défaut pour s'assurer que f_provider est initialisé
                f_provider = (x_val) -> -uExact.getSecondDerivative(x_val);
                break;
        }

        // Construction du système Ax = d_rhs
        switch (type) {
            case TYPE1: // -u'' + u = f  => (-1/h^2)u_i-1 + (2/h^2 + 1)u_i + (-1/h^2)u_i+1 = f_i
                for (int i = 0; i < n-1; i++) {
                    int actual_idx = i + 1;
                    a_sub[i] = -1.0/(h*h);
                    b_diag[i] = 2.0/(h*h) + 1.0;
                    c_sur[i] = -1.0/(h*h);
                    d_rhs[i] = f_provider.apply(x[actual_idx]); // x_coords -> x
                }
                // Ajustement pour les conditions aux limites u0 et u1
                d_rhs[0] -= a_sub[0] * u0;
                if (n-1 > 0) a_sub[0] = 0;

                if (n-1 > 0) d_rhs[n-2] -= c_sur[n-2] * u1;
                if (n-2 >=0 && n-2 < c_sur.length) c_sur[n-2] = 0;
                break;
                
            case TYPE2: // -u'' + u' = f. (-u_i-1 + 2u_i - u_i+1)/h^2 + (u_i+1 - u_i-1)/(2h) = f_i
                        // (-1/h^2 - 1/(2h))u_i-1 + (2/h^2)u_i + (-1/h^2 + 1/(2h))u_i+1 = f_i
                for (int i = 0; i < n-1; i++) {
                    int actual_idx = i + 1;
                    a_sub[i] = -1.0/(h*h) - 1.0/(2*h);
                    b_diag[i] = 2.0/(h*h);
                    c_sur[i] = -1.0/(h*h) + 1.0/(2*h);
                    d_rhs[i] = f_provider.apply(x[actual_idx]); // x_coords -> x
                }
                d_rhs[0] -= a_sub[0] * u0;
                if (n-1 > 0) a_sub[0] = 0;

                if (n-1 > 0) d_rhs[n-2] -= c_sur[n-2] * u1;
                if (n-2 >=0 && n-2 < c_sur.length) c_sur[n-2] = 0;
                break;
                
            case TYPE3: // -u'' = f  =>  -u_i-1 + 2u_i - u_i+1 = h^2 * f_i
                for (int i = 0; i < n-1; i++) {
                    int actual_idx = i + 1;
                    a_sub[i] = -1.0;
                    b_diag[i] = 2.0;
                    c_sur[i] = -1.0;
                    d_rhs[i] = h*h * f_provider.apply(x[actual_idx]); // x_coords -> x
                }
                d_rhs[0] += u0;
                if (n-1 > 0) a_sub[0] = 0;
                if (n-1 > 0) d_rhs[n-2] += u1;
                if (n-2 >=0 && n-2 < c_sur.length) c_sur[n-2] = 0;
                break;
        }
        
        // Résolution du système tridiagonal pour les points intérieurs u_inner = [u_1, ..., u_{n-1}]
        double[] u_inner = solveTridiagonal(a_sub, b_diag, c_sur, d_rhs);
        
        // Construction de la solution numérique complète (avec conditions aux limites u0, u1)
        double[] numerical = new double[n+1];
        double[] analytical = new double[n+1]; // Sera remplie par les valeurs de u(x) exact
        
        numerical[0] = u0;
        numerical[n] = u1;
        for (int i = 0; i < u_inner.length; i++) { // u_inner a n-1 éléments
            numerical[i+1] = u_inner[i];
        }
        
        // Solution analytique (u(x) exacte)
        for (int i = 0; i <= n; i++) {
            // La fonction getAnalyticalSolutionValue retourne la solution analytique u(x)
            // pour le problème posé avec les CLs u0, u1 et f(x) dérivé de uExact.
            analytical[i] = getAnalyticalSolutionValue(x[i], uExact, type, u0, u1);
        }
        
        // Calcul de l'erreur L-infini
        double max_abs_error = 0;
        for (int i = 0; i <= n; i++) {
            double current_abs_error = Math.abs(numerical[i] - analytical[i]);
            if (current_abs_error > max_abs_error) {
                max_abs_error = current_abs_error;
            }
        }
        // La variable 'error' dans la classe Solution stockera maintenant l'erreur L-infini
        
        return new Solution(x, numerical, analytical, max_abs_error, n, type, uExact); // x_coords -> x
    }
    
    // Calcul de l'ordre de convergence
    public static double calculateConvergenceOrder(double error1, double error2, int n1, int n2) {
        return Math.log(error1 / error2) / Math.log((double)n2 / n1);
    }
    
    // Classe pour la visualisation graphique
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
            
            // Axes
            g2.setColor(Color.BLACK);
            g2.drawLine(margin, height + margin, width + margin, height + margin); // x-axis
            g2.drawLine(margin, margin, margin, height + margin); // y-axis
            
            if (showError) {
                drawErrorGraph(g2, margin, width, height);
            } else {
                drawSolutionGraph(g2, margin, width, height);
            }
        }
        
        private void drawSolutionGraph(Graphics2D g2, int margin, int width, int height) {
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
            
            // Légende
            g2.setColor(Color.BLACK);
            g2.setFont(new Font("Arial", Font.BOLD, 16));
            // Affiche le type d'équation et la solution exacte u(x) utilisée pour dériver f(x)
            g2.drawString(sol.type.getDescription() + ", u_{exact}: " + sol.exactSolution.getName(), margin, 25);
            g2.drawString(String.format("N = %d, CL: u(0)=%.1f, u(1)=%.1f", sol.n, sol.numerical[0], sol.numerical[sol.n]), margin, 45);
            
            g2.setFont(new Font("Arial", Font.PLAIN, 12));
            g2.setColor(Color.BLUE);
            g2.drawString("— Solution analytique", width - 150, margin + 20);
            g2.setColor(Color.RED);
            g2.drawString("--- Solution numérique", width - 150, margin + 40);
            g2.setColor(Color.BLACK);
            g2.drawString(String.format("Erreur L∞: %.2e", sol.error), width - 150, margin + 60);
        }
        
        private void drawErrorGraph(Graphics2D g2, int margin, int width, int height) {
            if (solutions.size() < 2) return;
            
            // Filtrer les erreurs nulles ou négatives avant le log
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
            
            // Titre
            g2.setColor(Color.BLACK);
            g2.setFont(new Font("Arial", Font.BOLD, 16));
            Solution firstSol = solutions.get(0);
            g2.drawString("Évolution de l'erreur L∞ - " + firstSol.type.getDescription(), margin, 25);
            g2.drawString("Solution exacte u(x): " + firstSol.exactSolution.getName(), margin, 45); // Modifié pour afficher le nom de u(x)
            
            // Labels des axes
            g2.setFont(new Font("Arial", Font.PLAIN, 12));
            g2.drawString("Nombre de mailles N (log)", width/2, height + margin + 40); // Changé "mailles" en "mailles N"
            
            // Rotation pour le label Y
            AffineTransform orig = g2.getTransform();
            g2.rotate(-Math.PI/2);
            g2.drawString("Erreur L∞ (log)", -height/2 - 50, 20); // Changé L2 en L∞
            g2.setTransform(orig);
        }
    }
    
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            // Choix de la solution exacte u(x) et du type d'équation
            UValueProvider[] exactSolutions = {new USinPiX(), new UXCube()}; // Modifié
            EquationType[] types = {EquationType.TYPE3, EquationType.TYPE1, EquationType.TYPE2}; // TYPE3 en premier
            
            JFrame mainFrame = new JFrame("Résolveur d'équations différentielles 1D (Différences Finies)"); // Titre mis à jour
            mainFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            mainFrame.setLayout(new BorderLayout());
            
            JPanel controlPanel = new JPanel();
            JComboBox<UValueProvider> exactSolutionCombo = new JComboBox<>(exactSolutions); // Modifié
            JComboBox<EquationType> typeCombo = new JComboBox<>(types);
            JButton solveButton = new JButton("Résoudre et Afficher Solution");
            JButton errorButton = new JButton("Afficher Courbe d'Erreur");
            
            // Renderer pour UValueProvider dans JComboBox
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
             // Renderer pour EquationType dans JComboBox (déjà bon si toString est bien défini, sinon similaire)
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

            controlPanel.add(new JLabel("Solution exacte u(x):")); // Modifié
            controlPanel.add(exactSolutionCombo);
            controlPanel.add(new JLabel("Type d'équation:")); // Modifié
            controlPanel.add(typeCombo);
            controlPanel.add(solveButton);
            controlPanel.add(errorButton);
            
            JTabbedPane tabbedPane = new JTabbedPane();
            mainFrame.add(controlPanel, BorderLayout.NORTH);
            mainFrame.add(tabbedPane, BorderLayout.CENTER);
            
            // Conditions aux limites
            final double u0 = 0.0;
            final double u1 = 1.0;

            solveButton.addActionListener(e -> {
                UValueProvider selectedExactSolution = (UValueProvider) exactSolutionCombo.getSelectedItem();
                EquationType selectedType = (EquationType) typeCombo.getSelectedItem();
                
                tabbedPane.removeAll();
                
                int[] meshSizesForDisplay = {10, 20, 40, 80}; // Pour affichage individuel, moins de tabs
                System.out.println("\n=== Solutions pour " + selectedType.getDescription() + " avec " + selectedExactSolution.getName() + " ===");
                System.out.println("Conditions aux limites: u(0)=" + u0 + ", u(1)=" + u1);

                for (int n : meshSizesForDisplay) {
                    Solution sol = solve(n, selectedType, selectedExactSolution, u0, u1);
                    // Affichage graphique
                    GraphPanel panel = new GraphPanel(Arrays.asList(sol), false);
                    tabbedPane.addTab("Solution (N = " + n + ")", panel);
                    // Log console
                    System.out.printf("N = %d: Erreur L∞ = %.6e\n", n, sol.error);
                }
            });

            errorButton.addActionListener(e -> {
                UValueProvider selectedExactSolution = (UValueProvider) exactSolutionCombo.getSelectedItem();
                EquationType selectedType = (EquationType) typeCombo.getSelectedItem();

                tabbedPane.removeAll(); // Optionnel: ou ouvrir une nouvelle fenêtre pour la courbe d'erreur

                // Résolution pour différents nombres de mailles pour la courbe d'erreur
                int[] meshSizesForErrorCurve = {10, 20, 40, 80, 160, 320}; // Tailles demandées
                List<Solution> solutions = new ArrayList<>();
                
                System.out.println("\n=== Calcul de la courbe d'erreur pour " + selectedType.getDescription() + " avec " + selectedExactSolution.getName() + " ===");
                System.out.println("Conditions aux limites: u(0)=" + u0 + ", u(1)=" + u1);
                
                for (int n : meshSizesForErrorCurve) {
                    Solution sol = solve(n, selectedType, selectedExactSolution, u0, u1);
                    solutions.add(sol);
                    
                    System.out.printf("N = %d: Erreur L∞ = %.6e\n", n, sol.error);
                    
                    // Calcul de l'ordre de convergence
                    if (solutions.size() > 1) {
                        Solution prevSol = solutions.get(solutions.size() - 2);
                        double order = calculateConvergenceOrder(prevSol.error, sol.error, prevSol.n, sol.n);
                        System.out.printf("Ordre de convergence: %.2f\n", order);
                    }
                    
                    // Ajouter un onglet pour chaque solution
                    GraphPanel panel = new GraphPanel(Arrays.asList(sol), false);
                    tabbedPane.addTab("n = " + n, panel);
                }
                
                System.out.println();
                
                // Affichage de la courbe d'erreur dans un nouvel onglet ou une nouvelle fenêtre
                GraphPanel errorPanel = new GraphPanel(solutions, true); // 'solutions' contient déjà les résultats pour meshSizesForErrorCurve
                // Option 1: Ajouter comme onglet
                // tabbedPane.addTab("Courbe d'Erreur L∞", errorPanel);
                // Option 2: Nouvelle fenêtre (comme c'était avant, mais avec les bonnes données)
                JFrame errorFrame = new JFrame("Évolution de l'erreur L∞");
                errorFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE); // Important pour ne pas quitter toute l'appli
                errorFrame.add(errorPanel);
                errorFrame.pack(); // Ajuste la taille au contenu
                errorFrame.setLocationRelativeTo(mainFrame);
                errorFrame.setVisible(true);
            });
            
            // Suppression du deuxième ActionListener redondant pour errorButton qui était ici.

            mainFrame.pack(); // Ajuste la taille de la fenêtre principale aux composants
            mainFrame.setLocationRelativeTo(null);
            mainFrame.setVisible(true);
            
            // Sélection initiale et déclenchement pour exemple
            exactSolutionCombo.setSelectedIndex(0); // USinPiX
            typeCombo.setSelectedIndex(0);        // TYPE3
            // Déclencher l'affichage de la courbe d'erreur initialement peut-être ?
            // Ou laisser l'utilisateur cliquer. Pour l'instant, pas de clic auto.
            // solveButton.doClick(); // Affiche les solutions individuelles
            // errorButton.doClick(); // Affiche la courbe d'erreur
        });
    }
}