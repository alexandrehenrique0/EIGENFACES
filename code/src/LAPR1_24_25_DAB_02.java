import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class LAPR1_24_25_DAB_02 {
    public static final int MAX_SIZE_ROWS = 256;
    public static final int MAX_SIZE_COLS = 256;
    public static final int MIN_SIZE_ROWS = 1;
    public static final int MIN_SIZE_COLS = 1;
    public static Scanner SCANNER = new Scanner(System.in);
    public static Scanner SCANNER_CSV = new Scanner(System.in);
    public static Scanner SCANNER_IMAGE = new Scanner(System.in);

    public static void main(String[] args) {
        int function = 0;
        int vectorNumbers = 0;
        String csvLocation = "";
        String imageFolderLocation = "";
        if (check_Correct_Parameters(args)) {
            function = receive_Function(args);
            vectorNumbers = receive_Number_Vectors(args);
            csvLocation = receive_CSV_Location(args);
            imageFolderLocation = receive_Image_Location(args);

            try {
                SCANNER_CSV = new Scanner(new File(csvLocation));
                File[] fileList = new File(imageFolderLocation).listFiles();
            } catch (FileNotFoundException e) {
                System.out.println("Erro ao abrir os arquivos: " + e.getMessage());
                System.exit(1);
            }
            double[][] matrixCSVDouble = get_Matrix_From_CSV(csvLocation);
        } else if (args.length == 0) {
            // Mostrar as opções num menu e receber os parâmetros
            ui_Function_Parameter_Menu();
            function = receive_Function(null);

            ui_Vector_Numbers_Parameter_Menu();
            vectorNumbers = receive_Number_Vectors(null);

            ui_CSV_Location_Parameter_Menu();
            csvLocation = receive_CSV_Location(null);

            ui_Image_Location_Parameter_Menu();
            imageFolderLocation = receive_Image_Location(null);
            // ---------------------------------------

            try {
                SCANNER_CSV = new Scanner(new File(csvLocation));
                File[] fileList = new File(imageFolderLocation).listFiles();
            } catch (FileNotFoundException e) {
                System.out.println("Erro ao abrir os arquivos: " + e.getMessage());
                System.exit(1);
            }
            double[][] matrixCSVDouble = get_Matrix_From_CSV(csvLocation);

        } else {
            System.out.println("Parâmetros inválidos");
            System.exit(1);
        }
    }

    //* Verificações de parâmetros
    public static boolean check_Correct_Parameters(String[] parameters) {
        if (parameters.length == 8) {
            return parameters[0].equals("-f") && parameters[2].equals("-k") && parameters[4].equals("-i") && parameters[6].equals("-j");
        }
        return false;
    }
    public static boolean check_function(int function) {
        return function >= 1 && function <= 3;
    }
    public static boolean check_csvLocation(String csvLocation) {
        File csv = new File(csvLocation);
        if (csvLocation.equals("")) {
            return false;
        } else if (!csvLocation.contains(".csv")) {
            return false;
        } else return csv.exists();
    }
    public static boolean check_imageFolderLocation(String imageFolderLocation) {
        File imageDirectory = new File(imageFolderLocation);
        if (imageFolderLocation.isEmpty()) {
            return false;
        }
        return imageDirectory.exists();
    }
    //* ---------------------------------------

    // Menus de opções
    public static void ui_Function_Parameter_Menu() {
        System.out.println("------------- Que função deseja realizar? -------------");
        System.out.println("1 - Decomposição Própria de uma Matriz Simétrica");
        System.out.println("2 - Reconstrução de Imagens usando Eigenfaces");
        System.out.println("3 - Identificação de imagem mais próxima");
        System.out.println("-------------------------------------------------------");
        System.out.printf("Opção: ");
    }

    public static void ui_Vector_Numbers_Parameter_Menu() {
        System.out.println("------ Quantos vetores próprios deseja utilizar? ------");
        System.out.printf("Quantidade: ");
    }

    public static void ui_CSV_Location_Parameter_Menu() {
        System.out.println("---- Qual a localização do csv que deseja utilizar? ---");
        System.out.printf("Localização: ");
    }

    public static void ui_Image_Location_Parameter_Menu() {
        System.out.println("-------- Qual a localização da base de imagens? -------");
        System.out.printf("Localização: ");
    }

    // Receber parâmetros
    public static int receive_Function(String[] args) {
        int functionArgs = 0;
        if (args == null) {
            int function = SCANNER.nextInt();
            if (!check_function(function)) {
                error_Invalid_Option();
            }
            return function;
        } else {
            functionArgs = Integer.parseInt(args[1]);
            if (!check_function(functionArgs)) {
                error_Invalid_Option();
            }
            return functionArgs;
        }
    }
    public static int receive_Number_Vectors(String[] args) {
        int vectorNumbersArgs = 0;
        if (args == null) {
            int vectorNumbers = SCANNER.nextInt();
            return vectorNumbers;
        } else {
            vectorNumbersArgs = Integer.parseInt(args[3]);
            return vectorNumbersArgs;
        }
    }
    public static String receive_CSV_Location(String[] args) {
        String csvLocationArgs = "";
        if (args == null) {
            String csvLocation = SCANNER.next();
            if (!check_csvLocation(csvLocation)) {
                error_Location_Not_Found();
            }
            return csvLocation;
        } else {
            csvLocationArgs = args[5];
            if (!check_csvLocation(csvLocationArgs)) {
                error_Location_Not_Found();
            }
            return csvLocationArgs;
        }
    }

    public static String receive_Image_Location(String[] args) {
        String imageFolderLocationArgs = "";
        if (args == null) {
            String imageFolderLocation = SCANNER.next();
            if (!check_imageFolderLocation(imageFolderLocation)) {
               error_Location_Not_Found();
            }
            return imageFolderLocation;
        } else {
            imageFolderLocationArgs = args[7];
            if (!check_imageFolderLocation(imageFolderLocationArgs)) {
                error_Location_Not_Found();
            }
            return imageFolderLocationArgs;
        }
    }

    // ---------------------------------------

    //* Leitura de CSV
    public static double[][] get_Matrix_From_CSV(String csvLocation) {
        int[] dimensions = get_Dimensions();
        int rows = dimensions[0];
        int cols = dimensions[1];

        double[][] matrix = new double[rows][cols];
        populate_Matrix(matrix, csvLocation);

        return matrix;
    }
    private static int[] get_Dimensions() {
        int rows = 0;
        int cols = 0;
        while (SCANNER_CSV.hasNextLine()) {
            String line = SCANNER_CSV.nextLine().trim();
            if (!line.isEmpty()) {
                if (rows == 0) {
                    cols = line.split(",").length;
                }
                rows++;
            }
        }
        if (check_Size_Boundaries(rows, cols)) {
            throw new IllegalArgumentException("Erro: A matriz ultrapassa a dimensão máxima permitida de 256x256. Dimensão atual: " + rows + "x" + cols);
        }
        return new int[]{rows, cols};
    }

    private static boolean check_Size_Boundaries(int rows, int cols) {
        return rows > MAX_SIZE_ROWS || cols > MAX_SIZE_COLS || rows < MIN_SIZE_ROWS || cols < MIN_SIZE_COLS;
    }

    private static void populate_Matrix(double[][] matrix, String csvLocation) {
        try {
            SCANNER_CSV = new Scanner(new File(csvLocation));
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Erro ao reabrir o arquivo CSV: " + e.getMessage());
        }
        int row = 0;
        while (SCANNER_CSV.hasNextLine()) {
            String line = SCANNER_CSV.nextLine().trim();
            if (!line.isEmpty()) {
                populate_Row(matrix, row, line);
                row++;
            }
        }
        SCANNER_CSV.close();
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

    //* ---------------------------------------


    //* Decomposição de matriz simétrica
    public static EigenDecomposition decomposeMatrix(double[][] arrayParaDecompor) {
        Array2DRowRealMatrix matrix = new Array2DRowRealMatrix(arrayParaDecompor);
        EigenDecomposition eigenDecomposition = new EigenDecomposition(matrix);

        // demonstração de como utilizar os resultados da decomposição
//        RealMatrix eigenVectors = eigenDecomposition.getV();
//        RealMatrix eigenValues = eigenDecomposition.getD();
//        RealMatrix eigenVectorsTranspose = eigenDecomposition.getVT();

        return eigenDecomposition;
    }

    //* Metodo para printar matrizes corretamente no console.
    public static void printMatrix(RealMatrix matrix, String matrixName) {
        double[][] matrixToPrint = matrix.getData();
        System.out.println("Matriz: " + matrixName + " ↓");

        // Linha superior
        for (int i = 0; i < matrixToPrint[0].length; i++) {
            System.out.print("____________");
        }

        System.out.println(); // pula uma linha (estética)

        for (int j = 0; j < matrixToPrint.length; j++) {
            for (int i = 0; i < matrixToPrint[0].length; i++) {
                if (i == 0) {
                    System.out.printf("|%8.3f\t", matrixToPrint[j][i]);
                } else if (i == matrixToPrint[0].length - 1) {
                    System.out.printf("%8.3f\t|", matrixToPrint[j][i]);
                } else {
                    System.out.printf("%8.3f\t", matrixToPrint[j][i]);
                }
            }
            System.out.println();
        }
        // Linha final
        for (int i = 0; i < matrixToPrint[0].length; i++) {
            System.out.print("============");
        }
        System.out.println(); // pula duas linha (estética)
        System.out.println();
    }

    //! Error Messages
    public static void error_Location_Not_Found() {
        System.out.println("Localização inválida ou ficheiro não existe");
        System.exit(1);
    }
    public static void error_Invalid_Option() {
        System.out.println("Opção inválida");
        System.exit(1);
    }
    //! Error Messages
}