package finite_difference_gauss_seidel_serial;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class SerialGaussSeidel extends JFrame {

    private static final int CELL_SIZE = 100;
    private static final int GRID_OFFSET = 50;
    private int N = 3; // Grid size (NxN)
    private double[][] u; // Solution array
    private double h = 1.0 / (N - 1); // Step size, assuming domain [0,1]x[0,1]
    private double[][] exactSolution;
    private double[][] errorMatrix;
    private List<Double> l2ErrorHistory; // For L2 error vs iteration

    public SerialGaussSeidel() {
        super("Serial Gauss-Seidel 3x3");
        u = new double[N][N];
        exactSolution = new double[N][N];
        errorMatrix = new double[N][N];
        l2ErrorHistory = new ArrayList<>();

        initializeBoundaryConditions();
        calculateExactSolution();

        // GUI Setup
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLayout(new BorderLayout(10, 10));

        // Solve first, then update GUI
        solve(100, 1e-6); // Default solve parameters. For 3x3, converges very fast.
        calculateErrorMatrix();

        JPanel heatmapContainer = new JPanel(new GridLayout(1, 3, 10, 10));
        heatmapContainer.add(new HeatmapPanel(u, "Numerical Solution U"));
        heatmapContainer.add(new HeatmapPanel(exactSolution, "Exact Solution 2x+y"));
        heatmapContainer.add(new HeatmapPanel(errorMatrix, "Absolute Error |U - Exact|"));

        add(heatmapContainer, BorderLayout.CENTER);

        LinePlotPanel convergencePlot = new LinePlotPanel(l2ErrorHistory, "L2 Error vs. Iteration");
        add(convergencePlot, BorderLayout.SOUTH);

        // Adjust size - make it taller to accommodate the plot
        int heatmapWidth = 3 * (N * CELL_SIZE + 2 * GRID_OFFSET) + 100;
        int heatmapHeight = N * CELL_SIZE + 2 * GRID_OFFSET + 70; // Extra for titles etc.
        int plotHeight = 250; // Estimated height for the plot panel
        setSize(heatmapWidth, heatmapHeight + plotHeight);
        setLocationRelativeTo(null); // Center the window
    }

    private void calculateExactSolution() {
        for (int row = 0; row < N; row++) {
            for (int col = 0; col < N; col++) {
                double y_val = row * h;
                double x_val = col * h;
                exactSolution[row][col] = 2 * x_val + y_val;
            }
        }
    }

    // Inner class for drawing a line plot
    private static class LinePlotPanel extends JPanel {
        private final List<Double> dataPoints;
        private final String plotTitle;
        private static final int PADDING = 30;
        private static final int POINT_RADIUS = 3;

        public LinePlotPanel(List<Double> dataPoints, String plotTitle) {
            this.dataPoints = (dataPoints == null) ? new ArrayList<>() : dataPoints;
            this.plotTitle = plotTitle;
            setBorder(BorderFactory.createTitledBorder(plotTitle));
            setPreferredSize(new Dimension(400, 200)); // Default size
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g;
            g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            if (dataPoints.isEmpty()) {
                g2d.drawString("No data to plot.", getWidth() / 2 - 40, getHeight() / 2);
                return;
            }

            int width = getWidth();
            int height = getHeight();

            // Find min/max of data for y-axis scaling
            double minVal = dataPoints.stream().min(Double::compare).orElse(0.0);
            double maxVal = dataPoints.stream().max(Double::compare).orElse(1.0);
            if (minVal == maxVal) { // Adjust if all values are the same
                minVal -= 0.5 * Math.abs(minVal); // Add some padding
                maxVal += 0.5 * Math.abs(maxVal);
                 if (minVal == maxVal) { // If still same (e.g. all zeros)
                    minVal = -1.0; maxVal = 1.0;
                 }
            }


            // X-axis (Iterations)
            g2d.drawLine(PADDING, height - PADDING, width - PADDING, height - PADDING);
            // Y-axis (Value)
            g2d.drawLine(PADDING, PADDING, PADDING, height - PADDING);

            // Labels for axes
            FontMetrics fm = g2d.getFontMetrics();
            String xLabel = "Iteration";
            g2d.drawString(xLabel, width / 2 - fm.stringWidth(xLabel) / 2, height - PADDING / 2 + fm.getAscent()/2);

            // Y-axis label requires rotation
            String yLabel = "L2 Error Norm"; // Or a generic "Value"
            int yLabelWidth = fm.stringWidth(yLabel);
            g2d.rotate(-Math.PI / 2);
            g2d.drawString(yLabel, -height/2 - yLabelWidth/2 , PADDING / 2 + fm.getAscent()/2 -5 );
            g2d.rotate(Math.PI / 2);


            // Plot points and lines
            g2d.setColor(Color.BLUE.darker());
            double xScale = (dataPoints.size() <= 1) ? (width - 2 * PADDING) : (double) (width - 2 * PADDING) / (dataPoints.size() - 1);
            double yScale = (maxVal - minVal == 0) ? (height - 2 * PADDING) : (double) (height - 2 * PADDING) / (maxVal - minVal);

            int lastX = -1, lastY = -1;

            for (int i = 0; i < dataPoints.size(); i++) {
                int xVal = PADDING + (int) (i * xScale);
                int yVal = height - PADDING - (int) ((dataPoints.get(i) - minVal) * yScale);

                g2d.fillOval(xVal - POINT_RADIUS, yVal - POINT_RADIUS, 2 * POINT_RADIUS, 2 * POINT_RADIUS);
                if (lastX != -1) {
                    g2d.drawLine(lastX, lastY, xVal, yVal);
                }
                lastX = xVal;
                lastY = yVal;
            }

            // Y-axis scale labels
            g2d.setColor(Color.BLACK);
            g2d.drawString(String.format("%.2e", maxVal), PADDING - fm.stringWidth(String.format("%.2e", maxVal)) - 5, PADDING + fm.getAscent()/2);
            g2d.drawString(String.format("%.2e", minVal), PADDING - fm.stringWidth(String.format("%.2e", minVal)) - 5, height - PADDING + fm.getAscent()/2);

            // X-axis scale labels
            g2d.drawString("0", PADDING - fm.stringWidth("0")/2, height - PADDING + fm.getAscent() + 3);
            if (dataPoints.size() > 1) {
                String endIterLabel = String.valueOf(dataPoints.size() - 1);
                g2d.drawString(endIterLabel, width - PADDING - fm.stringWidth(endIterLabel)/2, height - PADDING + fm.getAscent() + 3);
            }
        }
    }

    private void calculateErrorMatrix() {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                errorMatrix[i][j] = Math.abs(u[i][j] - exactSolution[i][j]);
            }
        }
    }


    private void initializeBoundaryConditions() {
        // u[row][col] = 2*x + y where x = col*h, y = row*h
        for (int row = 0; row < N; row++) {
            for (int col = 0; col < N; col++) {
                if (row == 0 || row == N - 1 || col == 0 || col == N - 1) { // Boundary points
                    double y_val = row * h;
                    double x_val = col * h;
                    u[row][col] = 2 * x_val + y_val;
                } else {
                    u[row][col] = 0; // Initial guess for interior points
                }
            }
        }
    }

    private double calculateL2NormError(double[][] approxSolution, double[][] exactSol) {
        double sumSqError = 0.0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                // Consider only interior points for error calculation if desired,
                // or all points. For now, all points.
                double error = approxSolution[i][j] - exactSol[i][j];
                sumSqError += error * error;
            }
        }
        // The L2 norm is sqrt(h^2 * sum_sq_error).
        // For simplicity here, or if h=1 effectively for error comparison across iterations,
        // we can use sqrt(sum_sq_error / (N*N)) (RMS error) or just sqrt(sum_sq_error).
        // Let's use sqrt(sum_sq_error) as a measure of total error magnitude.
        // A more formal L2 norm would be Math.sqrt(sumSqError * h * h);
        return Math.sqrt(sumSqError / (N * N)); // RMS error
    }

    public void solve(int maxIterations, double tolerance) {
        l2ErrorHistory.clear();
        // u[row][col]
        for (int k = 0; k < maxIterations; k++) {
            double maxDiff = 0.0;
            // Iterate over interior points: row from 1 to N-2, col from 1 to N-2
            for (int row = 1; row < N - 1; row++) {
                for (int col = 1; col < N - 1; col++) {
                    double old_u_ij = u[row][col];
                    // u_ij = (u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1}) / 4
                    // u[row][col] = (u[row+1][col] + u[row-1][col] + u[row][col+1] + u[row][col-1]) / 4.0
                    u[row][col] = (u[row + 1][col] + u[row - 1][col] + u[row][col + 1] + u[row][col - 1]) / 4.0;
                    maxDiff = Math.max(maxDiff, Math.abs(u[row][col] - old_u_ij));
                }
            }

            l2ErrorHistory.add(calculateL2NormError(u, exactSolution));

            if (maxDiff < tolerance) {
                System.out.println("Converged after " + (k + 1) + " iterations.");
                return;
            }
        }
        System.out.println("Reached maximum iterations (" + maxIterations + ").");
    }

    // Inner class for drawing heatmap
    private static class HeatmapPanel extends JPanel {
        private final double[][] data;
        private final String title;
        private double globalMin = Double.MAX_VALUE;
        private double globalMax = Double.MIN_VALUE;

        public HeatmapPanel(double[][] data, String title) {
            this.data = data;
            this.title = title;
            setBorder(BorderFactory.createTitledBorder(title));

            for (double[] row : data) {
                for (double val : row) {
                    if (val < globalMin) globalMin = val;
                    if (val > globalMax) globalMax = val;
                }
            }
            // Handle case where all values are the same
            if (globalMin == globalMax) {
                globalMax = globalMin + 1.0; // Avoid division by zero in getColor
            }
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2d = (Graphics2D) g;
            int panelWidth = getWidth();
            int panelHeight = getHeight();
            int gridSize = data.length;

            // Calculate cell size based on panel dimensions to fit the title and some padding
            int titleHeight = 30; // Approximate height for titled border
            int actualCellSize = Math.min((panelWidth - GRID_OFFSET) / gridSize, (panelHeight - GRID_OFFSET - titleHeight) / gridSize);


            for (int row = 0; row < gridSize; row++) {
                for (int col = 0; col < gridSize; col++) {
                    double value = data[row][col];
                    float normalized = (globalMax - globalMin == 0) ? 0.5f : (float) ((value - globalMin) / (globalMax - globalMin));
                    // Blue (low) to Red (high) colormap
                    Color color = new Color(normalized, 0, 1 - normalized);
                    g2d.setColor(color);
                    g2d.fillRect(GRID_OFFSET / 2 + col * actualCellSize, GRID_OFFSET / 2 + titleHeight + row * actualCellSize, actualCellSize, actualCellSize);
                    g2d.setColor(Color.BLACK);
                    g2d.drawRect(GRID_OFFSET / 2 + col * actualCellSize, GRID_OFFSET / 2 + titleHeight + row * actualCellSize, actualCellSize, actualCellSize);

                    // Draw text
                    String text = String.format("%.2f", value);
                    FontMetrics fm = g2d.getFontMetrics();
                    Rectangle2D textBounds = fm.getStringBounds(text, g2d);
                    int textX = GRID_OFFSET / 2 + col * actualCellSize + (actualCellSize - (int) textBounds.getWidth()) / 2;
                    int textY = GRID_OFFSET / 2 + titleHeight + row * actualCellSize + (actualCellSize - (int) textBounds.getHeight()) / 2 + fm.getAscent();

                    // Determine text color based on background brightness
                    double luminance = 0.299 * color.getRed() + 0.587 * color.getGreen() + 0.114 * color.getBlue();
                    if (luminance < 128) {
                        g2d.setColor(Color.WHITE);
                    } else {
                        g2d.setColor(Color.BLACK);
                    }
                    g2d.drawString(text, textX, textY);
                }
            }
        }

        @Override
        public Dimension getPreferredSize() {
            return new Dimension(data.length * CELL_SIZE + GRID_OFFSET, data.length * CELL_SIZE + GRID_OFFSET + 30);
        }
    }


    public void printSolution() {
        System.out.println("Solution u:");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                System.out.printf("%.4f\t", u[i][j]);
            }
            System.out.println();
        }
    }

    public double[][] getSolution() {
        return u;
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            SerialGaussSeidel solverFrame = new SerialGaussSeidel();
            solverFrame.setVisible(true);
        });
    }
}
