
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class Rafael {
    public static Scanner SCANNER;
    static {
        try {
            SCANNER = new Scanner(new File("input.txt"));
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }
    }

    public static double[][] get_Matrix_From_CSV(){
        int[] dimensions = get_Dimensions();
        int rows = dimensions[0];
        int cols = dimensions[1];

        double[][] matrix = new double[rows][cols];
        populate_Matrix(matrix);

        return matrix;
    }
    private static int[] get_Dimensions(){
        int rows = 0;
        int cols = 0;
        while (SCANNER.hasNextLine()) {
            String line = SCANNER.nextLine().trim();
            if (!line.isEmpty()) {
                if (rows == 0) {
                    cols = line.split(",").length;
                }
                rows++;
            }
        }
        SCANNER.close();
        return new int[]{rows, cols};
    }
    private static void populate_Matrix(double[][] matrix) {
        int row = 0;
        while (SCANNER.hasNextLine()) {
            String line = SCANNER.nextLine().trim();
            if (!line.isEmpty()) {
                populate_Row(matrix, row, line);
                row++;
            }
        }
        SCANNER.close();
    }

    private static void populate_Row(double[][] matrix, int row, String line) {
        String[] values = line.split(",");
        for (int col = 0; col < values.length; col++) {
            try {
                matrix[row][col] = Double.parseDouble(values[col].trim());
            } catch (NumberFormatException e) {
                matrix[row][col] = 0; // or any default value you prefer
            }
        }
    }
}