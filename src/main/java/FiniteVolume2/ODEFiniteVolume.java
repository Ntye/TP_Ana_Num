package FiniteVolume2;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Classe principale pour la résolution d'équations différentielles par volumes finis
 * Résout -u'' = f avec f = sin(π*x) + y^3
 */
public class ODEFiniteVolume extends JFrame {
    private static final double PI = Math.PI;
    private static final int WINDOW_SIZE = 800;
    private int currentN = 10;
    private double[][] currentSolution;
    private double[][] exactSolution;
    private double[] meshSizes = {10, 20, 40};
    private double[] errors = new double[meshSizes.length];
    private JPanel mainPanel;
    private JLabel infoLabel;
    
    public ODEFiniteVolume() {
        setTitle("Résolution EDP par Volumes Finis - Méthodes Comparatives");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLayout(new BorderLayout());
        
        mainPanel = new JPanel();
        mainPanel.setLayout(new GridLayout(2, 2));
        
        infoLabel = new JLabel("Cliquez sur les graphiques pour voir les valeurs");
        infoLabel.setHorizontalAlignment(SwingConstants.CENTER);
        infoLabel.setFont(new Font("Arial", Font.BOLD, 14));
        
        add(infoLabel, BorderLayout.NORTH);
        add(mainPanel, BorderLayout.CENTER);
        
        // Panel de contrôle
        JPanel controlPanel = new JPanel();
        JButton btn10 = new JButton("N=10");
        JButton btn20 = new JButton("N=20");
        JButton btn40 = new JButton("N=40");
        JButton btnError = new JButton("Courbe d'erreur");
        
        btn10.addActionListener(e -> updateDisplay(10));
        btn20.addActionListener(e -> updateDisplay(20));
        btn40.addActionListener(e -> updateDisplay(40));
        btnError.addActionListener(e -> showErrorCurve());
        
        controlPanel.add(btn10);
        controlPanel.add(btn20);
        controlPanel.add(btn40);
        controlPanel.add(btnError);
        
        add(controlPanel, BorderLayout.SOUTH);
        
        // Calcul initial
        runComparativeStudy();
        updateDisplay(10);
        
        setSize(WINDOW_SIZE, WINDOW_SIZE + 100);
        setLocationRelativeTo(null);
    }
    
    /**
     * Fonction source f(x,y) = sin(π*x) + y^3
     */
    private double sourceFunction(double x, double y) {
        return Math.sin(PI * x) + Math.pow(y, 3);
    }
    
    /**
     * Solution exacte approximée (pour validation)
     */
    private double exactSolutionFunction(double x, double y) {
        // Solution approximée basée sur la série de Fourier
        return Math.sin(PI * x) / (PI * PI) + y * y * y * y * y / 20.0;
    }
    
    /**
     * Résolution par méthode directe (Elimination de Gauss)
     */
    public class DirectSolver {
        public double[][] solve(int n) {
            double h = 1.0 / (n + 1);
            int size = n * n;
            double[][] A = new double[size][size];
            double[] b = new double[size];
            
            // Construction de la matrice A et du vecteur b
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    int idx = i * n + j;
                    double x = (j + 1) * h;
                    double y = (i + 1) * h;
                    
                    // Coefficient diagonal
                    A[idx][idx] = 4.0 / (h * h);
                    
                    // Coefficients adjacents
                    if (i > 0) A[idx][idx - n] = -1.0 / (h * h);
                    if (i < n - 1) A[idx][idx + n] = -1.0 / (h * h);
                    if (j > 0) A[idx][idx - 1] = -1.0 / (h * h);
                    if (j < n - 1) A[idx][idx + 1] = -1.0 / (h * h);
                    
                    b[idx] = sourceFunction(x, y);
                }
            }
            
            // Résolution par élimination de Gauss
            double[] solution = gaussElimination(A, b);
            
            // Conversion en matrice 2D
            double[][] result = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    result[i][j] = solution[i * n + j];
                }
            }
            
            return result;
        }
        
        private double[] gaussElimination(double[][] A, double[] b) {
            int n = A.length;
            
            // Elimination avant
            for (int k = 0; k < n - 1; k++) {
                for (int i = k + 1; i < n; i++) {
                    if (Math.abs(A[k][k]) < 1e-10) continue;
                    double factor = A[i][k] / A[k][k];
                    for (int j = k; j < n; j++) {
                        A[i][j] -= factor * A[k][j];
                    }
                    b[i] -= factor * b[k];
                }
            }
            
            // Substitution arrière
            double[] x = new double[n];
            for (int i = n - 1; i >= 0; i--) {
                x[i] = b[i];
                for (int j = i + 1; j < n; j++) {
                    x[i] -= A[i][j] * x[j];
                }
                if (Math.abs(A[i][i]) > 1e-10) {
                    x[i] /= A[i][i];
                }
            }
            
            return x;
        }
    }
    
    /**
     * Résolution par Gauss-Seidel avec relaxation
     */
    public class GaussSeidelSolver {
        private double omega = 1.5; // Facteur de relaxation
        private int maxIter = 1000;
        private double tolerance = 1e-6;
        
        public double[][] solve(int n) {
            double h = 1.0 / (n + 1);
            double[][] u = new double[n][n];
            double[][] uOld = new double[n][n];
            
            // Initialisation
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    u[i][j] = 0.0;
                }
            }
            
            // Itérations Gauss-Seidel avec relaxation
            for (int iter = 0; iter < maxIter; iter++) {
                // Sauvegarde solution précédente
                for (int i = 0; i < n; i++) {
                    System.arraycopy(u[i], 0, uOld[i], 0, n);
                }
                
                // Mise à jour des valeurs
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        double x = (j + 1) * h;
                        double y = (i + 1) * h;
                        
                        double sum = 0.0;
                        if (i > 0) sum += u[i-1][j];
                        if (i < n-1) sum += uOld[i+1][j];
                        if (j > 0) sum += u[i][j-1];
                        if (j < n-1) sum += uOld[i][j+1];
                        
                        double newVal = (h * h * sourceFunction(x, y) + sum) / 4.0;
                        u[i][j] = (1 - omega) * uOld[i][j] + omega * newVal;
                    }
                }
                
                // Test de convergence
                double maxDiff = 0.0;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        maxDiff = Math.max(maxDiff, Math.abs(u[i][j] - uOld[i][j]));
                    }
                }
                
                if (maxDiff < tolerance) {
                    System.out.println("Gauss-Seidel convergé en " + iter + " itérations");
                    break;
                }
            }
            
            return u;
        }
    }
    
    /**
     * Résolution par Gauss-Seidel parallélisé
     */
    public class ParallelGaussSeidelSolver {
        private double omega = 1.5;
        private int maxIter = 1000;
        private double tolerance = 1e-6;
        private ExecutorService executor = Executors.newFixedThreadPool(4);
        
        public double[][] solve(int n) {
            double h = 1.0 / (n + 1);
            double[][] u = new double[n][n];
            double[][] uOld = new double[n][n];
            
            // Initialisation
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    u[i][j] = 0.0;
                }
            }
            
            // Itérations parallèles
            for (int iter = 0; iter < maxIter; iter++) {
                // Sauvegarde solution précédente
                for (int i = 0; i < n; i++) {
                    System.arraycopy(u[i], 0, uOld[i], 0, n);
                }
                
                // Mise à jour parallèle par damier (red-black)
                updateRedBlack(u, uOld, n, h, true);  // Points rouges
                updateRedBlack(u, uOld, n, h, false); // Points noirs
                
                // Test de convergence
                double maxDiff = 0.0;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        maxDiff = Math.max(maxDiff, Math.abs(u[i][j] - uOld[i][j]));
                    }
                }
                
                if (maxDiff < tolerance) {
                    System.out.println("Gauss-Seidel parallèle convergé en " + iter + " itérations");
                    break;
                }
            }
            
            executor.shutdown();
            return u;
        }
        
        private void updateRedBlack(double[][] u, double[][] uOld, int n, double h, boolean red) {
            Future<?>[] futures = new Future[4];
            
            for (int t = 0; t < 4; t++) {
                final int thread = t;
                futures[t] = executor.submit(() -> {
                    int startRow = thread * n / 4;
                    int endRow = (thread + 1) * n / 4;
                    
                    for (int i = startRow; i < endRow; i++) {
                        for (int j = 0; j < n; j++) {
                            // Damier: rouge si (i+j) pair, noir sinon
                            if (((i + j) % 2 == 0) != red) continue;
                            
                            double x = (j + 1) * h;
                            double y = (i + 1) * h;
                            
                            double sum = 0.0;
                            if (i > 0) sum += u[i-1][j];
                            if (i < n-1) sum += uOld[i+1][j];
                            if (j > 0) sum += u[i][j-1];
                            if (j < n-1) sum += uOld[i][j+1];
                            
                            double newVal = (h * h * sourceFunction(x, y) + sum) / 4.0;
                            u[i][j] = (1 - omega) * uOld[i][j] + omega * newVal;
                        }
                    }
                });
            }
            
            // Attendre la fin de tous les threads
            for (Future<?> future : futures) {
                try {
                    future.get();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }
    
    /**
     * Calcul de l'erreur L2
     */
    private double calculateError(double[][] numerical, double[][] exact) {
        double error = 0.0;
        int n = numerical.length;
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double diff = numerical[i][j] - exact[i][j];
                error += diff * diff;
            }
        }
        
        return Math.sqrt(error / (n * n));
    }
    
    /**
     * Calcul de l'ordre de convergence
     */
    private double calculateConvergenceOrder(double[] errors) {
        if (errors.length < 2) return 0.0;
        
        double order = 0.0;
        for (int i = 1; i < errors.length; i++) {
            if (errors[i] > 0 && errors[i-1] > 0) {
                order += Math.log(errors[i-1] / errors[i]) / Math.log(2.0);
            }
        }
        
        return order / (errors.length - 1);
    }
    
    /**
     * Étude comparative des trois méthodes
     */
    private void runComparativeStudy() {
        DirectSolver directSolver = new DirectSolver();
        GaussSeidelSolver gsSolver = new GaussSeidelSolver();
        ParallelGaussSeidelSolver pgsSolver = new ParallelGaussSeidelSolver();
        
        System.out.println("=== ÉTUDE COMPARATIVE ===");
        System.out.println("Équation: -u'' = sin(π*x) + y³");
        System.out.println();
        
        for (int idx = 0; idx < meshSizes.length; idx++) {
            int n = (int) meshSizes[idx];
            System.out.println("Maillage " + n + "x" + n + ":");
            
            // Solution exacte
            double[][] exact = new double[n][n];
            double h = 1.0 / (n + 1);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double x = (j + 1) * h;
                    double y = (i + 1) * h;
                    exact[i][j] = exactSolutionFunction(x, y);
                }
            }
            
            // Méthode directe
            long startTime = System.currentTimeMillis();
            double[][] directSol = directSolver.solve(n);
            long directTime = System.currentTimeMillis() - startTime;
            double directError = calculateError(directSol, exact);
            
            // Gauss-Seidel
            startTime = System.currentTimeMillis();
            double[][] gsSol = gsSolver.solve(n);
            long gsTime = System.currentTimeMillis() - startTime;
            double gsError = calculateError(gsSol, exact);
            
            // Gauss-Seidel parallèle
            startTime = System.currentTimeMillis();
            double[][] pgsSol = pgsSolver.solve(n);
            long pgsTime = System.currentTimeMillis() - startTime;
            double pgsError = calculateError(pgsSol, exact);
            
            errors[idx] = directError; // Utilisation de l'erreur directe pour la courbe
            
            System.out.printf("  Méthode directe: Erreur=%.6f, Temps=%dms%n", directError, directTime);
            System.out.printf("  Gauss-Seidel: Erreur=%.6f, Temps=%dms%n", gsError, gsTime);
            System.out.printf("  GS Parallèle: Erreur=%.6f, Temps=%dms%n", pgsError, pgsTime);
            System.out.println();
        }
        
        double convergenceOrder = calculateConvergenceOrder(errors);
        System.out.printf("Ordre de convergence: %.3f%n", convergenceOrder);
        System.out.println("==================");
    }
    
    /**
     * Mise à jour de l'affichage
     */
    private void updateDisplay(int n) {
        currentN = n;
        
        // Calcul des solutions
        DirectSolver solver = new DirectSolver();
        currentSolution = solver.solve(n);
        
        exactSolution = new double[n][n];
        double h = 1.0 / (n + 1);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double x = (j + 1) * h;
                double y = (i + 1) * h;
                exactSolution[i][j] = exactSolutionFunction(x, y);
            }
        }
        
        // Mise à jour des panneaux
        mainPanel.removeAll();
        
        mainPanel.add(new SolutionPanel("Solution Numérique", currentSolution));
        mainPanel.add(new SolutionPanel("Solution Exacte", exactSolution));
        mainPanel.add(new ErrorPanel("Erreur Absolue", currentSolution, exactSolution));
        mainPanel.add(new InfoPanel());
        
        mainPanel.revalidate();
        mainPanel.repaint();
        
        double error = calculateError(currentSolution, exactSolution);
        infoLabel.setText(String.format("N=%d, Erreur L2=%.6f", n, error));
    }
    
    /**
     * Affichage de la courbe d'erreur
     */
    private void showErrorCurve() {
        JFrame errorFrame = new JFrame("Évolution de l'erreur");
        errorFrame.setSize(600, 400);
        errorFrame.setLocationRelativeTo(this);
        
        JPanel errorPanel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                Graphics2D g2d = (Graphics2D) g;
                g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
                
                int width = getWidth() - 80;
                int height = getHeight() - 80;
                int offsetX = 40;
                int offsetY = 40;
                
                // Axes
                g2d.drawLine(offsetX, offsetY + height, offsetX + width, offsetY + height);
                g2d.drawLine(offsetX, offsetY, offsetX, offsetY + height);
                
                // Étiquettes
                g2d.drawString("Nombre de mailles", offsetX + width/2, offsetY + height + 30);
                g2d.rotate(-Math.PI/2);
                g2d.drawString("Erreur L2", -offsetY - height/2, 20);
                g2d.rotate(Math.PI/2);
                
                // Courbe d'erreur
                g2d.setColor(Color.RED);
                g2d.setStroke(new BasicStroke(2));
                
                for (int i = 0; i < errors.length - 1; i++) {
                    int x1 = offsetX + (int)(i * width / (errors.length - 1));
                    int x2 = offsetX + (int)((i + 1) * width / (errors.length - 1));
                    
                    double maxError = 0.0;
                    for (double error : errors) {
                        maxError = Math.max(maxError, error);
                    }
                    
                    int y1 = offsetY + height - (int)(errors[i] * height / maxError);
                    int y2 = offsetY + height - (int)(errors[i + 1] * height / maxError);
                    
                    g2d.drawLine(x1, y1, x2, y2);
                    g2d.fillOval(x1 - 3, y1 - 3, 6, 6);
                }
                
                // Dernier point
                int lastX = offsetX + width;
                int lastY = offsetY + height - (int)(errors[errors.length - 1] * height / getMaxError());
                g2d.fillOval(lastX - 3, lastY - 3, 6, 6);
                
                // Valeurs
                g2d.setColor(Color.BLACK);
                for (int i = 0; i < errors.length; i++) {
                    int x = offsetX + (int)(i * width / (errors.length - 1));
                    g2d.drawString(String.valueOf((int)meshSizes[i]), x - 10, offsetY + height + 15);
                    g2d.drawString(String.format("%.2e", errors[i]), x - 20, offsetY + height - (int)(errors[i] * height / getMaxError()) - 10);
                }
            }
            
            private double getMaxError() {
                double max = 0.0;
                for (double error : errors) {
                    max = Math.max(max, error);
                }
                return max;
            }
        };
        
        errorFrame.add(errorPanel);
        errorFrame.setVisible(true);
    }
    
    /**
     * Panel pour afficher une solution avec carte de couleur
     */
    private class SolutionPanel extends JPanel implements MouseListener, MouseMotionListener {
        private String title;
        private double[][] data;
        private double minVal, maxVal;
        
        public SolutionPanel(String title, double[][] data) {
            this.title = title;
            this.data = data;
            this.setBorder(BorderFactory.createTitledBorder(title));
            this.addMouseListener(this);
            this.addMouseMotionListener(this);
            
            // Calcul min/max
            minVal = maxVal = data[0][0];
            for (double[] row : data) {
                for (double val : row) {
                    minVal = Math.min(minVal, val);
                    maxVal = Math.max(maxVal, val);
                }
            }
        }
        
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            
            int width = getWidth() - 20;
            int height = getHeight() - 40;
            int n = data.length;
            
            double cellWidth = (double) width / n;
            double cellHeight = (double) height / n;
            
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double normalizedVal = (data[i][j] - minVal) / (maxVal - minVal);
                    Color color = getColor(normalizedVal);
                    
                    g.setColor(color);
                    int x = 10 + (int)(j * cellWidth);
                    int y = 30 + (int)(i * cellHeight);
                    g.fillRect(x, y, (int)cellWidth + 1, (int)cellHeight + 1);
                }
            }
            
            // Grille
            g.setColor(Color.GRAY);
            for (int i = 0; i <= n; i++) {
                int x = 10 + (int)(i * cellWidth);
                int y = 30 + (int)(i * cellHeight);
                g.drawLine(x, 30, x, 30 + height);
                g.drawLine(10, y, 10 + width, y);
            }
            
            // Échelle de couleur
            drawColorScale(g);
        }
        
        private Color getColor(double value) {
            // Palette de couleur bleue vers rouge
            if (value < 0.5) {
                return new Color(0, (int)(255 * 2 * value), 255);
            } else {
                return new Color((int)(255 * 2 * (value - 0.5)), 255, (int)(255 * 2 * (1 - value)));
            }
        }
        
        private void drawColorScale(Graphics g) {
            int scaleWidth = 20;
            int scaleHeight = 100;
            int x = getWidth() - scaleWidth - 10;
            int y = getHeight() - scaleHeight - 10;
            
            for (int i = 0; i < scaleHeight; i++) {
                double val = (double) i / scaleHeight;
                g.setColor(getColor(val));
                g.drawLine(x, y + scaleHeight - i, x + scaleWidth, y + scaleHeight - i);
            }
            
            g.setColor(Color.BLACK);
            g.drawRect(x, y, scaleWidth, scaleHeight);
            g.drawString(String.format("%.3f", maxVal), x + scaleWidth + 5, y + 5);
            g.drawString(String.format("%.3f", minVal), x + scaleWidth + 5, y + scaleHeight);
        }
        
        @Override
        public void mouseClicked(MouseEvent e) {
            showValueAtPoint(e.getX(), e.getY());
        }
        
        @Override
        public void mouseMoved(MouseEvent e) {
            showValueAtPoint(e.getX(), e.getY());
        }
        
        private void showValueAtPoint(int mouseX, int mouseY) {
            int width = getWidth() - 20;
            int height = getHeight() - 40;
            int n = data.length;
            
            double cellWidth = (double) width / n;
            double cellHeight = (double) height / n;
            
            int i = (int)((mouseY - 30) / cellHeight);
            int j = (int)((mouseX - 10) / cellWidth);
            
            if (i >= 0 && i < n && j >= 0 && j < n) {
                double h = 1.0 / (n + 1);
                double x = (j + 1) * h;
                double y = (i + 1) * h;
                
                String info = String.format("%s - Point(%.3f,%.3f): %.6f", 
                    title, x, y, data[i][j]);
                infoLabel.setText(info);
            }
        }
        
        @Override public void mousePressed(MouseEvent e) {}
        @Override public void mouseReleased(MouseEvent e) {}
        @Override public void mouseEntered(MouseEvent e) {}
        @Override public void mouseExited(MouseEvent e) {}
        @Override public void mouseDragged(MouseEvent e) {}
    }
    
    /**
     * Panel pour afficher l'erreur
     */
    private class ErrorPanel extends SolutionPanel {
        public ErrorPanel(String title, double[][] numerical, double[][] exact) {
            super(title, calculateErrorMatrix(numerical, exact));
        }
        
        private static double[][] calculateErrorMatrix(double[][] numerical, double[][] exact) {
            int n = numerical.length;
            double[][] error = new double[n][n];
            
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    error[i][j] = Math.abs(numerical[i][j] - exact[i][j]);
                }
            }
            
            return error;
        }
    }
    
    /**
     * Panel d'informations
     */
    private class InfoPanel extends JPanel {
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            
            g.setColor(Color.BLACK);
            g.setFont(new Font("Arial", Font.BOLD, 12));
            
            int y = 30;
            g.drawString("Méthodes de résolution:", 10, y);
            y += 20;
            g.drawString("• Méthode directe (Élimination de Gauss)", 10, y);
            y += 15;
            g.drawString("• Gauss-Seidel avec relaxation", 10, y);
            y += 15;
            g.drawString("• Gauss-Seidel parallélisé", 10, y);
            y += 25;
            
            g.drawString("Équation: -u'' = sin(π*x) + y³", 10, y);
            y += 20;
            
            if (currentSolution != null && exactSolution != null) {
                double error = calculateError(currentSolution, exactSolution);
                g.drawString(String.format("Erreur L2: %.6f", error), 10, y);
                y += 15;
                
                g.drawString("Ordre de convergence: " + 
                    String.format("%.3f", calculateConvergenceOrder(errors)), 10, y);
            }
        }
    }
    
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            try {
                UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
            } catch (Exception e) {
                e.printStackTrace();
            }
            
            new ODEFiniteVolume().setVisible(true);
        });
    }
}