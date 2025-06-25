package finite_difference_gauss_seidel_parallel;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveAction;
import java.util.concurrent.atomic.DoubleAccumulator;

public class ParallelGaussSeidel extends JFrame {

    private static final int N_GRID_SIZE = 100; // Grid size (NxN)
    private static final int MAX_CELL_SIZE_FOR_TEXT = 20; // Max cell size in pixels to attempt drawing text
    private double[][] u; // Solution array
    private double h = 1.0 / (N_GRID_SIZE - 1); // Step size, assuming domain [0,1]x[0,1]
    private ForkJoinPool forkJoinPool;
    private double[][] exactSolution;
    private double[][] errorMatrix;
    private List<Double> maxDiffHistory = new ArrayList<>(); // To store maxDiff per iteration

    // GUI Components
    private HeatmapPanel numericalSolutionPanel;
    private HeatmapPanel exactSolutionPanel;
    private HeatmapPanel errorPanel;
    // private JPanel convergencePlotPanel; // Will be added in the next step
    private ConvergencePlotPanel convergencePlotPanelComponent;


    public ParallelGaussSeidel() {
        super("Parallel Gauss-Seidel " + N_GRID_SIZE + "x" + N_GRID_SIZE);
        u = new double[N_GRID_SIZE][N_GRID_SIZE];
        exactSolution = new double[N_GRID_SIZE][N_GRID_SIZE];
        errorMatrix = new double[N_GRID_SIZE][N_GRID_SIZE];

        initializeBoundaryConditions();
        calculateExactSolution();

        // Use a ForkJoinPool with a number of threads equal to available processors
        int numProcessors = Runtime.getRuntime().availableProcessors();
        System.out.println("Using " + numProcessors + " processors for ForkJoinPool.");
        forkJoinPool = new ForkJoinPool(numProcessors);

        // GUI Setup
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        // Layout will be set later, potentially with a plot panel
        // For now, let's plan for 3 heatmaps + 1 plot area
        // Main panel for heatmaps
        JPanel heatmapContainer = new JPanel(new GridLayout(1, 3, 5, 5));

        // Solve first, then update GUI
        System.out.println("Solving... This might take a moment for N=" + N_GRID_SIZE);
        long startTime = System.currentTimeMillis();
        solve(20000, 1e-7); // Default solve parameters from previous main
        long endTime = System.currentTimeMillis();
        System.out.println("Solved in " + (endTime - startTime) + " ms.");

        calculateErrorMatrix();

        heatmapContainer.add(new HeatmapPanel(u, "Numerical Solution U"));
        heatmapContainer.add(new HeatmapPanel(exactSolution, "Exact Solution (2x+y)"));
        heatmapContainer.add(new HeatmapPanel(errorMatrix, "Absolute Error |U - Exact|"));

        // Placeholder for convergence plot panel
        JPanel plotPanel = new JPanel(); // Will be replaced by actual plot
        plotPanel.setBorder(BorderFactory.createTitledBorder("Convergence (maxDiff vs Iteration)"));


        // Main content pane
        JPanel mainContentPanel = new JPanel(new BorderLayout());
        mainContentPanel.add(heatmapContainer, BorderLayout.CENTER);
        // mainContentPanel.add(plotPanel, BorderLayout.SOUTH); // Add plot later
        convergencePlotPanelComponent = new ConvergencePlotPanel(maxDiffHistory);
        mainContentPanel.add(convergencePlotPanelComponent, BorderLayout.SOUTH);


        numericalSolutionPanel = new HeatmapPanel(u, "Numerical Solution U");
        exactSolutionPanel = new HeatmapPanel(exactSolution, "Exact Solution (2x+y)");
        errorPanel = new HeatmapPanel(errorMatrix, "Absolute Error |U - Exact|");

        heatmapContainer.add(numericalSolutionPanel);
        heatmapContainer.add(exactSolutionPanel);
        heatmapContainer.add(errorPanel);

        setContentPane(mainContentPanel);
        pack(); // Adjusts window size to preferred sizes of components.
        // setSize(1200, 650); // Adjusted for plot, or rely on pack()
        setLocationRelativeTo(null); // Center the window
    }

    private void calculateExactSolution() {
        // u[row][col] = 2*x + y where x = col*h, y = row*h
        for (int row = 0; row < N_GRID_SIZE; row++) {
            for (int col = 0; col < N_GRID_SIZE; col++) {
                double y_val = row * h;
                double x_val = col * h;
                exactSolution[row][col] = 2 * x_val + y_val;
            }
        }
    }

    private static class ConvergencePlotPanel extends JPanel {
        private final List<Double> history;
        private static final int PADDING = 25;
        private static final int POINT_RADIUS = 2;

        public ConvergencePlotPanel(List<Double> history) {
            this.history = history;
            setBorder(BorderFactory.createTitledBorder("Convergence: log(maxDiff) vs. Iteration"));
            setPreferredSize(new Dimension(800, 200)); // Initial preferred size
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            if (history == null || history.isEmpty()) {
                g2d.drawString("No convergence data available.", getWidth() / 2 - 50, getHeight() / 2);
                return;
            }

            int width = getWidth();
            int height = getHeight();

            // Find min/max of log(history) for y-axis scaling
            double minLogDiff = Double.MAX_VALUE;
            double maxLogDiff = Double.MIN_VALUE;
            for (Double diff : history) {
                if (diff > 0) { // Avoid log(0) or log(negative)
                    double logDiff = Math.log10(diff);
                    if (logDiff < minLogDiff) minLogDiff = logDiff;
                    if (logDiff > maxLogDiff) maxLogDiff = logDiff;
                }
            }

            if (minLogDiff == Double.MAX_VALUE) { // All diffs were <=0 or list empty after filtering
                 g2d.drawString("Convergence data not plottable (all diffs <= 0).", getWidth() / 2 - 100, getHeight() / 2);
                return;
            }
            // If min and max are the same (e.g. only one point), add some padding for scale
             if (minLogDiff == maxLogDiff) {
                minLogDiff -= 1;
                maxLogDiff += 1;
            }


            // X-axis (Iterations)
            g2d.drawLine(PADDING, height - PADDING, width - PADDING, height - PADDING);
            // Y-axis (log(maxDiff))
            g2d.drawLine(PADDING, PADDING, PADDING, height - PADDING);

            // Labels
            g2d.drawString("Iterations", width / 2 - 20, height - PADDING / 2);
            g2d.rotate(-Math.PI / 2);
            g2d.drawString("log10(maxDiff)", -height / 2 - 30, PADDING / 2 + 5);
            g2d.rotate(Math.PI / 2);

            // Plot points
            g2d.setColor(Color.BLUE);
            double xScale = (double) (width - 2 * PADDING) / (history.size() -1 == 0 ? 1 : history.size() -1) ;
            double yScale = (double) (height - 2 * PADDING) / (maxLogDiff - minLogDiff);

            int lastX = -1, lastY = -1;

            for (int i = 0; i < history.size(); i++) {
                if (history.get(i) <= 0) continue; // Skip non-positive diffs for log scale

                double logVal = Math.log10(history.get(i));
                int xVal = PADDING + (int) (i * xScale);
                int yVal = height - PADDING - (int) ((logVal - minLogDiff) * yScale);

                g2d.fillOval(xVal - POINT_RADIUS, yVal - POINT_RADIUS, 2 * POINT_RADIUS, 2 * POINT_RADIUS);
                if (lastX != -1) {
                    g2d.drawLine(lastX, lastY, xVal, yVal);
                }
                lastX = xVal;
                lastY = yVal;
            }

            // Y-axis scale labels
            g2d.setColor(Color.BLACK);
            FontMetrics fm = g2d.getFontMetrics();
            String minLabel = String.format("%.1e", Math.pow(10,minLogDiff));
            String maxLabel = String.format("%.1e", Math.pow(10,maxLogDiff));
            g2d.drawString(minLabel, PADDING - fm.stringWidth(minLabel) - 5, height - PADDING + fm.getAscent()/2);
            g2d.drawString(maxLabel, PADDING - fm.stringWidth(maxLabel) - 5, PADDING + fm.getAscent()/2 );

            // X-axis scale labels
            g2d.drawString("0", PADDING - fm.stringWidth("0")/2, height - PADDING + fm.getAscent() + 2);
            String endIterLabel = String.valueOf(history.size()-1);
            g2d.drawString(endIterLabel, width - PADDING - fm.stringWidth(endIterLabel)/2, height - PADDING + fm.getAscent() + 2);
        }
    }

    private void calculateErrorMatrix() {
        for (int row = 0; row < N_GRID_SIZE; row++) {
            for (int col = 0; col < N_GRID_SIZE; col++) {
                errorMatrix[row][col] = Math.abs(u[row][col] - exactSolution[row][col]);
            }
        }
    }

    private void initializeBoundaryConditions() {
        // u[row][col] = 2*x + y where x = col*h, y = row*h
        for (int row = 0; row < N_GRID_SIZE; row++) {
            for (int col = 0; col < N_GRID_SIZE; col++) {
                if (row == 0 || row == N_GRID_SIZE - 1 || col == 0 || col == N_GRID_SIZE - 1) { // Boundary points
                    double y_val = row * h;
                    double x_val = col * h;
                    u[row][col] = 2 * x_val + y_val;
                } else {
                    u[row][col] = 0.0; // Initial guess for interior points
                }
            }
        }
    }

    // Inner class for parallel computation task
    private class GaussSeidelTask extends RecursiveAction {
        private final int startRow, endRow, color; // 0 for red, 1 for black
        private final DoubleAccumulator maxDiffAccumulator;

        GaussSeidelTask(int startRow, int endRow, int color, DoubleAccumulator maxDiffAccumulator) {
            this.startRow = startRow;
            this.endRow = endRow;
            this.color = color;
            this.maxDiffAccumulator = maxDiffAccumulator;
        }

        @Override
        protected void compute() {
            // Threshold for sequential processing, ensuring startRow and endRow are valid interior indices
            int actualStartRow = Math.max(1, startRow);
            int actualEndRow = Math.min(N_GRID_SIZE - 1, endRow);

            if (actualEndRow - actualStartRow < 10) { // Adjust threshold as needed
                for (int i = actualStartRow; i < actualEndRow; i++) { // Iterate up to N_GRID_SIZE-2 for rows
                    for (int j = 1; j < N_GRID_SIZE - 1; j++) { // Iterate up to N_GRID_SIZE-2 for columns
                        if ((i + j) % 2 == color) { // Process only cells of the current color
                            double old_u_ij = u[i][j];
                            // u[row][col] = (u[row+1][col] + u[row-1][col] + u[row][col+1] + u[row][col-1]) / 4.0
                            u[i][j] = (u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1]) / 4.0;
                            maxDiffAccumulator.accumulate(Math.abs(u[i][j] - old_u_ij));
                        }
                    }
                }
            } else {
                int midRow = (actualStartRow + actualEndRow) / 2;
                invokeAll(new GaussSeidelTask(actualStartRow, midRow, color, maxDiffAccumulator),
                          new GaussSeidelTask(midRow, actualEndRow, color, maxDiffAccumulator));
            }
        }
    }

    public void solve(int maxIterations, double tolerance) {
        maxDiffHistory.clear(); // Clear history for a new solve run
        for (int k = 0; k < maxIterations; k++) {
            DoubleAccumulator maxDiffAccumulator = new DoubleAccumulator(Math::max, 0.0);

            // Red pass: iterate rows from 1 to N_GRID_SIZE-2 (inclusive for start, exclusive for end)
            GaussSeidelTask redTask = new GaussSeidelTask(1, N_GRID_SIZE - 1, 0, maxDiffAccumulator);
            forkJoinPool.invoke(redTask);

            // Black pass
            GaussSeidelTask blackTask = new GaussSeidelTask(1, N_GRID_SIZE - 1, 1, maxDiffAccumulator);
            forkJoinPool.invoke(blackTask);

            double maxDiff = maxDiffAccumulator.get();

            if (k % 10 == 0 || maxDiff < tolerance) { // Print progress less frequently for larger grids
                 System.out.printf("Iteration %d, Max Diff: %.6e\n", k + 1, maxDiff);
            }

            if (maxDiff < tolerance) {
                System.out.println("Converged after " + (k + 1) + " iterations.");
                forkJoinPool.shutdown();
                return;
            }
        }
        System.out.println("Reached maximum iterations (" + maxIterations + ").");
        forkJoinPool.shutdown();
    }

    public double[][] getSolution() {
        return u;
    }

    public List<Double> getMaxDiffHistory() {
        return maxDiffHistory;
    }

    // HeatmapPanel class (can be an inner class or a separate file)
    // For simplicity, defining it as an inner class here.
    // Adapted from SerialGaussSeidel's HeatmapPanel
    private static class HeatmapPanel extends JPanel {
        private final double[][] data;
        private final String title;
        private double dataMin = Double.MAX_VALUE;
        private double dataMax = Double.MIN_VALUE;
        private static final int PADDING = 20;
        private static final int LEGEND_WIDTH = 60;
        private static final int LEGEND_TEXT_PADDING = 5;


        public HeatmapPanel(double[][] data, String title) {
            this.data = data;
            this.title = title;
            this.setBorder(BorderFactory.createTitledBorder(title));
            calculateMinMax();
        }

        private void calculateMinMax() {
            dataMin = Double.MAX_VALUE;
            dataMax = Double.MIN_VALUE;
            for (double[] row : data) {
                for (double val : row) {
                    if (val < dataMin) dataMin = val;
                    if (val > dataMax) dataMax = val;
                }
            }
            if (dataMin == dataMax) { // Handle uniform data
                dataMax = dataMin + 1.0;
            }
        }

        public void updateData(double[][] newData) {
            // This method could be used if the panel needs to be updated dynamically
            // For now, data is set at construction
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g;

            int panelWidth = getWidth();
            int panelHeight = getHeight();
            int gridSize = data.length;

            int heatmapWidth = panelWidth - 2 * PADDING - LEGEND_WIDTH;
            int heatmapHeight = panelHeight - 2 * PADDING;

            double cellWidth = (double) heatmapWidth / gridSize;
            double cellHeight = (double) heatmapHeight / gridSize;

            for (int row = 0; row < gridSize; row++) {
                for (int col = 0; col < gridSize; col++) {
                    double value = data[row][col];
                    float normalized = (dataMax - dataMin == 0) ? 0.5f : (float) ((value - dataMin) / (dataMax - dataMin));
                    Color color = getColor(normalized);

                    g2d.setColor(color);
                    g2d.fillRect(PADDING + (int)(col * cellWidth), PADDING + (int)(row * cellHeight), (int)Math.ceil(cellWidth), (int)Math.ceil(cellHeight));
                }
            }
            drawColorScale(g2d, panelWidth - PADDING - LEGEND_WIDTH + LEGEND_TEXT_PADDING, PADDING, heatmapHeight);
        }

        private Color getColor(float value) { // value is normalized 0 to 1
            // Simple blue (low) to red (high) colormap
             return Color.getHSBColor( (1.0f - value) * 240f / 360f, 1.0f, 1.0f); // Hue from Blue to Red
        }

        private void drawColorScale(Graphics2D g2d, int x, int y, int height) {
            int barWidth = 20; // Width of the color bar
            for (int i = 0; i < height; i++) {
                float normalized = (float) i / (height -1); // from bottom (0) to top (1)
                g2d.setColor(getColor(1.0f - normalized)); // getColor expects 0=low, 1=high. Here, top is high.
                g2d.drawLine(x, y + i, x + barWidth, y + i);
            }
            g2d.setColor(Color.BLACK);
            g2d.drawRect(x, y, barWidth, height);

            g2d.drawString(String.format("%.2e", dataMax), x + barWidth + LEGEND_TEXT_PADDING, y + g2d.getFontMetrics().getAscent());
            g2d.drawString(String.format("%.2e", dataMin), x + barWidth + LEGEND_TEXT_PADDING, y + height);
             // Add mid-point value
            double midVal = dataMin + (dataMax - dataMin) / 2.0;
            g2d.drawString(String.format("%.2e", midVal), x + barWidth + LEGEND_TEXT_PADDING, y + height/2 + g2d.getFontMetrics().getAscent()/2);
        }

        @Override
        public Dimension getPreferredSize() {
            // Suggest a size, e.g., 400x400 for heatmap area + padding + legend
            int preferredHeatmapAreaSize = 300;
            return new Dimension(preferredHeatmapAreaSize + 2 * PADDING + LEGEND_WIDTH + 50, preferredHeatmapAreaSize + 2 * PADDING);
        }
    }


    public static void main(String[] args) {
         SwingUtilities.invokeLater(() -> {
            ParallelGaussSeidel solverFrame = new ParallelGaussSeidel();
            solverFrame.setVisible(true);

            // Optional: If ForkJoinPool is not shutdown by solver, shutdown on window close
            solverFrame.addWindowListener(new java.awt.event.WindowAdapter() {
                @Override
                public void windowClosing(java.awt.event.WindowEvent windowEvent) {
                    if (solverFrame.forkJoinPool != null && !solverFrame.forkJoinPool.isShutdown()) {
                        System.out.println("Shutting down ForkJoinPool from window listener.");
                        solverFrame.forkJoinPool.shutdown();
                    }
                }
            });
        });
    }
}
