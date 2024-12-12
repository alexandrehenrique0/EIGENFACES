
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Scanner;

public class Rafael {

    public static double[][] create_submatrix_remove_col(double[][] matrix, int col) {
        double[][] submatrix = new double[matrix.length][matrix[0].length - 1];
        int sub_i = 0;
        for (double[] doubles : matrix) {
            int sub_j = 0;
            for (int j = 0; j < matrix[0].length; j++) {
                if (j == col - 1) continue;
                submatrix[sub_i][sub_j] = doubles[j];
                sub_j++;
            }
            sub_i++;
        }
        return submatrix;
    }


}