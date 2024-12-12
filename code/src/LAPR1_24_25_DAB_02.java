import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;

public class LAPR1_24_25_DAB_02 {


    // CONSTANTES PARA LIMITES (MAX E MIN) DE MATRIZES
    public static final int MAX_SIZE_ROWS = 256;
    public static final int MAX_SIZE_COLS = 256;
    public static final int MIN_SIZE_ROWS = 1;
    public static final int MIN_SIZE_COLS = 1;

    // SCANNERS GLOBAIS
    public static Scanner SCANNER = new Scanner(System.in);
    public static Scanner SCANNER_CSV;
    public static Scanner SCANNER_IMAGE;
    public static File[] IMAGE_FILES;

    public static void main(String[] args) {

        // PARÂMETROS DE ENTRADA
        int function = 0;
        int vectorNumbers = 0;
        String csvLocation = "";
        String imageFolderLocation = "";

        // VERIFICAÇÃO DE PARÂMETROS (ARGS) PARA MÉTODOS DE INPUT
        if (check_Correct_Parameters_Structure(args)) {

            // Receber os parâmetros
            function = receive_Function(args);
            vectorNumbers = receive_Number_Vectors(args);
            csvLocation = receive_CSV_Location(args);
            imageFolderLocation = receive_Image_Location(args);

            // Verificar se os arquivo e diretório existem
            check_Existance_File_Directory(csvLocation, imageFolderLocation);

            // Obter a matriz do CSV
            double[][] matrixCSVDouble = get_Matrix_From_CSV(csvLocation);

            //FUnção que contém as funções principais
            switch_Primary_Functions(function, matrixCSVDouble, vectorNumbers);

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

            // Verificar se os arquivo e diretório existem
            check_Existance_File_Directory(csvLocation, imageFolderLocation);

            // Obter a matriz do CSV
            double[][] matrixCSVDouble = get_Matrix_From_CSV(csvLocation);

            //FUnção que contém as funções principais
            switch_Primary_Functions(function, matrixCSVDouble, vectorNumbers);

        } else {
            error_General("Erro: Parâmetros inválidos");
        }
    }

    public static void switch_Primary_Functions(int function, double[][] matrixCSVDouble, int vectorNumbers) {
        switch (function) {
            case 1:
                for (int i = 0; i < vectorNumbers; i++) {
                    //TODO Função Gabriel
                    //TODO double[][] matrixPk = create_subMatrix_remove_col(matriz P, col);

                }
                break;
            case 2:
                // TODO Reconstrução de Imagens usando Eigenfaces
                break;
            case 3:
                // TODO Identificação de imagem mais próxima
                break;
        }
    }

    //* ------------------ Verificações ------------------
    public static boolean check_Correct_Parameters_Structure(String[] parameters) {
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

    public static void check_Existance_File_Directory(String csvLocation, String imageFolderLocation) {
        try {
            SCANNER_CSV = new Scanner(new File(csvLocation));
            File[] IMAGE_FILES = new File(imageFolderLocation).listFiles();
        } catch (FileNotFoundException e) {
            error_General("Erro ao abrir os arquivos: " + e.getMessage());
        }
    }
    //* ------------------ Fim verificações ------------------


    //* ------------------ Menus de opções ------------------
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
    //* --------------------- Fim menus de opções ------------------


    //* ------------------ Receber parâmetros ------------------
    public static int receive_Function(String[] args) {
        int functionArgs = 0;
        if (args == null) {
            int function = SCANNER.nextInt();
            if (!check_function(function)) {
                error_General("Erro: Opção inválida");
            }
            return function;
        } else {
            functionArgs = Integer.parseInt(args[1]);
            if (!check_function(functionArgs)) {
                error_General("Erro: Opção inválida");
            }
            return functionArgs;
        }
    }

    public static int receive_Number_Vectors(String[] args) {
        int vectorNumbersArgs = 0;
        if (args == null) {
            return SCANNER.nextInt();
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
                error_General("Erro: Localização inválida csv");
            }
            return csvLocation;
        } else {
            csvLocationArgs = args[5];
            if (!check_csvLocation(csvLocationArgs)) {
                error_General("Erro: Localização inválida csv");
            }
            return csvLocationArgs;
        }
    }

    public static String receive_Image_Location(String[] args) {
        String imageFolderLocationArgs = "";
        if (args == null) {
            String imageFolderLocation = SCANNER.next();
            if (!check_imageFolderLocation(imageFolderLocation)) {
                error_General("Erro: Localização inválida da pasta de imagens");
            }
            return imageFolderLocation;
        } else {
            imageFolderLocationArgs = args[7];
            if (!check_imageFolderLocation(imageFolderLocationArgs)) {
                error_General("Erro: Localização inválida da pasta de imagens");
            }
            return imageFolderLocationArgs;
        }
    }
    //* -------------------- Fim receber parâmetros ------------------


    //* ------------------ Leitura de CSV ------------------
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
            error_General("Erro: Dimensões da matriz fora dos limites: " + rows + "x" + cols);
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
    //* ----------------------- Fim leitura de CSV ------------------


    //* ------------------ Operações com Matrizes ------------------
    public static double[][] sum_Matrices(double[][] matrizLeft, double[][] matrizRight) {
        double[][] matrizResultante = new double[matrizLeft.length][matrizLeft[0].length];
        for (int i = 0; i < matrizLeft.length; i++) {
            for (int j = 0; j < matrizLeft[0].length; j++) {
                matrizResultante[i][j] = matrizLeft[i][j] + matrizRight[i][j];
            }
        }
        return matrizResultante;
    }

    public static double[][] subtraction_Matrices(double[][] matrizLeft, double[][] matrizRight) {
        double[][] matrizResultante = new double[matrizLeft.length][matrizLeft[0].length];
        for (int i = 0; i < matrizLeft.length; i++) {
            for (int j = 0; j < matrizLeft[0].length; j++) {
                matrizResultante[i][j] = matrizLeft[i][j] - matrizRight[i][j];
            }
        }
        return matrizResultante;
    }

    public static double[][] multiply_Matrices(double[][] matrizLeft, double[][] matrizRight) {
        double[][] matrizResultante = new double[matrizLeft.length][matrizRight[0].length];
        for (int i = 0; i < matrizLeft.length; i++) {
            for (int j = 0; j < matrizRight[0].length; j++) {
                for (int k = 0; k < matrizRight.length; k++) {
                    matrizResultante[i][j] += matrizLeft[i][k] * matrizRight[k][j];
                }
            }
        }
        return matrizResultante;
    }

    public static double[][] multiply_Matrix_Escalar(double[][] matriz, double escalar) {
        double[][] matrizResultante = new double[matriz.length][matriz[0].length];
        for (int i = 0; i < matriz.length; i++) {
            for (int j = 0; j < matriz[0].length; j++) {
                matrizResultante[i][j] = matriz[i][j] * escalar;
            }
        }
        return matrizResultante;
    }

    public static double[][] transposed_Matrix(double[][] matriz) {
        double[][] matrizTransposta = new double[matriz[0].length][matriz.length];
        for (int i = 0; i < matriz.length; i++) {
            for (int j = 0; j < matriz[0].length; j++) {
                matrizTransposta[j][i] = matriz[i][j];
            }
        }
        return matrizTransposta;
    }

    //! Usar se preciso para remover apenas uma coluna
    public static double[][] create_subMatrix_remove_col(double[][] matrix, int col) {
        double[][] submatrix = new double[matrix.length][matrix[0].length - 1];
        int sub_i = 0;
        for (double[] doubles : matrix) {
            int sub_j = 0;
            for (int j = 0; j < matrix[0].length; j++) {
                if (j == col) continue;
                submatrix[sub_i][sub_j] = doubles[j];
                sub_j++;
            }
            sub_i++;
        }
        return submatrix;
    }

    //! -----------------------------------
    public static double[][] create_submatrix_Keep_cols(double[][] matrixP, double[][] keepColumns) {

        boolean[] keepColumnsBOL = new boolean[matrixP[0].length];

        for (double[] col : keepColumns) {
            keepColumnsBOL[(int) col[1]] = true;
        }

        double[][] submatrix = new double[matrixP.length][keepColumns.length];

        int sub_i = 0;
        for (double[] doubles : matrixP) {
            int sub_j = 0;
            for (int j = 0; j < doubles.length; j++) {
                if (keepColumnsBOL[j]) {
                    submatrix[sub_i][sub_j] = doubles[j];
                    sub_j++;
                }
            }
            sub_i++;
        }
        return submatrix;
    }
    //* ----------------- Fim operaçoes com matrizes ------------------


    //* -------------------- Printar Matrizes -----------------------
    // Metodo para printar matrizes corretamente no console.
    public static void print_Matrix(double[][] matrixToPrint, String matrixName) {
        System.out.println("Matriz: " + matrixName + " ↓");
        print_Line(matrixToPrint[0].length, "____________");

        for (double[] row : matrixToPrint) {
            System.out.print("|");
            for (int i = 0; i < row.length; i++) {
                System.out.printf("%8.3f\t", row[i]);
                if (i == row.length - 1) {
                    System.out.print("|");
                }
            }
            System.out.println();
        }
        print_Line(matrixToPrint[0].length, "============");
        System.out.println();
    }

    // Metodo para printar linhas de caracteres no console.
    private static void print_Line(int length, String pattern) {
        for (int i = 0; i < length; i++) {
            System.out.print(pattern);
        }
        System.out.println();
    }
    //* ----------------- Fim printar matrizes -----------------------


    //* -------------------- Funcionalidade 1 -----------------------
    // eu (gabriel) apenas colocarei os metodos que ainda não tem no main, depois podemos organizar melhor
    public static EigenDecomposition decompose_Matrix(double[][] arrayParaDecompor) {
        Array2DRowRealMatrix matrix = new Array2DRowRealMatrix(arrayParaDecompor);
        EigenDecomposition eigenDecomposition = new EigenDecomposition(matrix);

        /*? demonstração de como utilizar os resultados da decomposição
        RealMatrix eigenVectors = eigenDecomposition.getV();
        RealMatrix eigenValues = eigenDecomposition.getD();
        RealMatrix eigenVectorsTranspose = eigenDecomposition.getVT();
        ? ----------------------------------------------------------*/
        return eigenDecomposition;
    }

    public static double calculateEAM(double[][] A, double[][] Ak) {
        int M = A.length;
        int N = A[0].length;
        double erroAbsMed = 0;
        // Percorre cada elemento da matriz
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                erroAbsMed += Math.abs(A[i][j] - Ak[i][j]);
            }
        }

        // Calcula o erro médio
        return erroAbsMed / (M * N);
    }

    public static double[][] getEigenVectors(EigenDecomposition eigenDecomposition) {
        RealMatrix eigenVectors = eigenDecomposition.getV();
        return eigenVectors.getData();
    }

    public static double[][] getEigenValues(EigenDecomposition eigenDecomposition) {
        RealMatrix eigenValues = eigenDecomposition.getD();
        return eigenValues.getData();
    }

    public static double[][] getEigenVectorsTranspose(EigenDecomposition eigenDecomposition) {
        RealMatrix eigenVectorsTranspose = eigenDecomposition.getVT();
        return eigenVectorsTranspose.getData();
    }

    private static double[][] getEigenVectorsKArray(double[][] eigenVectorsArray, double[][] arrayValuesK) {
        double[][] newVectorsK = new double[eigenVectorsArray.length][arrayValuesK.length];
        //TODO ao inves de dar system exit, pode retornar a pergunta de quantos K o utilizador quer
        if (arrayValuesK.length > eigenVectorsArray[0].length) {
            System.out.println("O valor de k não pode ser maior que o número de colunas da matriz de autovetores.");
            System.exit(1);
        } else if (arrayValuesK.length < 1) {
            System.out.println("O valor de k não pode ser menor que 1.");
            System.exit(1);
        }

        for (int i = 0; i < eigenVectorsArray[0].length; i++) {
            for (int j = 0; j < arrayValuesK.length; j++) {
                newVectorsK[i][j] = eigenVectorsArray[i][(int) arrayValuesK[j][1]];
            }
        }


        return newVectorsK;
    }

    private static double[][] constructDiagonalMatrix(double[][] matrixvaluesK) {
        double[][] matrixvaluesKPrint = new double[matrixvaluesK.length][matrixvaluesK.length];
        for (int i = 0; i < matrixvaluesK.length; i++) {
            matrixvaluesKPrint[i][i] = matrixvaluesK[i][0];
        }
        return matrixvaluesKPrint;
    }

    private static double[][] getValuesAndIndexArray(double[][] eigenValuesArray, int k) {
        double[][] valuesAndIndexArray = new double[k][2];

        for (int i = 0; i < valuesAndIndexArray.length; i++) {
            valuesAndIndexArray[i][0] = Double.MIN_VALUE;
        }

        for (int i = 0; i < eigenValuesArray.length; i++) {
            double absValue = Math.abs(eigenValuesArray[i][i]);
            for (int j = 0; j < valuesAndIndexArray.length; j++) {
                if (absValue > Math.abs(valuesAndIndexArray[j][0])) {
                    for (int l = valuesAndIndexArray.length - 1; l > j; l--) {
                        valuesAndIndexArray[l][0] = valuesAndIndexArray[l - 1][0];
                        valuesAndIndexArray[l][1] = valuesAndIndexArray[l - 1][1];
                    }
                    valuesAndIndexArray[j][0] = eigenValuesArray[i][i];
                    valuesAndIndexArray[j][1] = i;
                    break;
                }
            }
        }

        return valuesAndIndexArray;
    }

    // serve para A e para Ak
    public static double[][] multiplyVectorsValuesVectorsTransposed(double[][] matrixVectors, double[][] matrixValues, double[][] matrixVectorsTranspose) {
        return multiply_Matrices(multiply_Matrices(matrixVectors, matrixValues), matrixVectorsTranspose);
    }

    //* ----------------- Fim funcionalidade 1 ------------------


    //! ------------------ Error Messages ------------------
    public static void error_General(String error) {
        System.out.println(error);
        System.exit(1);
    }

    //! ------------------ Fim error messages --------------
}