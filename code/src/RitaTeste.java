import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

    public class RitaTeste {

        // CONSTANTES PARA LIMITES (MAX E MIN) DE MATRIZES
        public static final int MAX_SIZE_ROWS = 256;
        public static final int MAX_SIZE_COLS = 256;
        public static final int MIN_SIZE_ROWS = 1;
        public static final int MIN_SIZE_COLS = 1;

        // SCANNERS GLOBAIS
        public static Scanner SCANNER = new Scanner(System.in);
        public static Scanner SCANNER_CSV;

        public static void main(String[] args) {

            // PARÂMETROS DE ENTRADA
            int function = 0;
            int vectorNumbers = 0;
            String csvLocation = "";
            String imageFolderLocation = "";

            boolean KeepRunning = true;

            do {


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

                    // Obter a matriz do CSV para a função 2
                    // TODO: metodo para ler matrizes e gerar uma matriz 3D
                    double[][][] matrixCSVDouble3D = get_Matrices_From_CSV_Folder(imageFolderLocation);

                    //FUnção que contém as funções principais
                    switch_Primary_Functions(function, matrixCSVDouble, vectorNumbers, csvLocation, matrixCSVDouble3D, imageFolderLocation);

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

                    // Obter a matriz do CSV para a função 2
                    // TODO: metodo para ler matrizes e gerar uma matriz 3D
                    double[][][] matrixCSVDouble3D = get_Matrices_From_CSV_Folder(imageFolderLocation);

                    //FUnção que contém as funções principais
                    switch_Primary_Functions(function, matrixCSVDouble, vectorNumbers, csvLocation, matrixCSVDouble3D, imageFolderLocation);


                } else {
                    error_General("Erro: Parâmetros inválidos");
                }
                System.out.println();
                System.out.println();
                System.out.println("-------------------------------------------------------------------------------------");
                System.out.print("Digite 'sair' para fechar o programa ou qualquer outra tecla para continuar: ");
                String input = SCANNER.next();
                System.out.println("-------------------------------------------------------------------------------------");
                if (input.equalsIgnoreCase("sair")) {
                    KeepRunning = false;
                }
            } while (KeepRunning);


        }

        public static void switch_Primary_Functions(int function, double[][] matrixCSVDouble, int vectorNumbers, String csvLocation, double[][][] matrixCSVDouble3D, String imageFolderLocation) {
            // Common variables
            String[] csvFiles = getCSVFileNames(imageFolderLocation);
            double[][] primaryMatrix = readCSVToArray(csvLocation);
            double[][][] imageMatrices = get_Matrices_From_CSV_Folder(imageFolderLocation);
            double[][] linearizedImages = new double[imageMatrices[0].length * imageMatrices[0].length][imageMatrices.length];
            double[] meanVector;
            double[][] phi;
            double[][] covariance;
            double[][] eigenVecs;
            double[][] selectedColumns;
            double[][] eigenVecsK;
            double[][] eigenfaces;
            double[][] normEigenfaces;
            double[][] weightsMatrix;
            double[] phiVector;
            double[] weightsVetorPrincipal;
            double[] euclideanDistances;
            int closestImageIndex;

            populateLinearizedImages(linearizedImages, imageMatrices);
            Scanner scanner = new Scanner(System.in);
            String entrada = ""; // Variável para armazenar a entrada do usuário

            // O loop continua enquanto a entrada não for "done"

            switch (function) {
                case 1:
                    print_Header_Function("Decomposição Própria de uma Matriz Simétrica");
                    EigenDecomposition eigenDecomposition = decompose_Matrix(matrixCSVDouble);
                    double[][] eigenVectors = getEigenVectors(eigenDecomposition);
                    double[][] eigenValues = getEigenValues(eigenDecomposition);
                    double[][] eigenVectorsTranspose = getEigenVectorsTranspose(eigenDecomposition);
                    double[][] valuesAndIndexArray = getValuesAndIndexArray(eigenValues, vectorNumbers);

                    double[][] newEigenVectorsK = create_submatrix_Keep_cols(eigenVectors, valuesAndIndexArray);
                    double[][] newEigenValuesK = constructDiagonalMatrix(valuesAndIndexArray);
                    double[][] newEigenVectorsTransposeK = transposed_Matrix(newEigenVectorsK);
                    double[][] resultingMatrixAk = multiplyVectorsValuesVectorsTransposed(newEigenVectorsK, newEigenValuesK, newEigenVectorsTransposeK);

                    double errorAbsMed = calculateMAE(matrixCSVDouble, resultingMatrixAk);
                    saveMatrixToFile(resultingMatrixAk, csvLocation, "Output/Func1");
                    print_Function_1(matrixCSVDouble, vectorNumbers, newEigenVectorsK, newEigenValuesK, newEigenVectorsTransposeK, resultingMatrixAk, errorAbsMed);
                    break;
                case 2:
                    print_Header_Function("Reconstrução de Imagens usando Eigenfaces");

                    meanVector = calculateMeanVector(linearizedImages);
                    phi = centralizeImages(linearizedImages, meanVector);
                    covariance = covariances(phi);
                    eigenVecs = eigenVectors(covariance);

                    if (vectorNumbers == -1 || vectorNumbers > eigenVecs[0].length) {
                        vectorNumbers = eigenVecs[0].length;
                    }

                    selectedColumns = getValuesAndIndexArray(eigenVecs, vectorNumbers);
                    eigenVecsK = create_submatrix_Keep_cols(eigenVecs, selectedColumns);
                    eigenfaces = multiply_Matrices(phi, eigenVecsK);
                    normEigenfaces = normalize(eigenfaces);

                    weightsMatrix = new double[phi.length][phi[0].length];
                    populateWeightsMatrix(weightsMatrix, phi, normEigenfaces);

                    for (int img = 0; img < linearizedImages[0].length; img++) {
                        double[] reconstructedImage = reconstructImage(meanVector, normEigenfaces, getColumn(weightsMatrix, img));
                        double[][] reconstructedImageMatrix = array1DToMatrix(reconstructedImage, imageMatrices[img]);

                        saveImage(reconstructedImageMatrix, csvFiles[img], "Output/Func2/ImagensReconstruidas");
                        saveMatrixToFile(reconstructedImageMatrix, csvFiles[img], "Output/Func2/Eigenfaces");
                    }
                    break;
                case 3:
                    print_Header_Function("Identificação de imagem mais próxima");

                    meanVector = calculateMeanVector(linearizedImages);
                    phi = centralizeImages(linearizedImages, meanVector);
                    covariance = covariances(phi);
                    eigenVecs = eigenVectors(covariance);

                    if (vectorNumbers == -1 || vectorNumbers > eigenVecs[0].length) {
                        vectorNumbers = eigenVecs[0].length;
                    }

                    selectedColumns = getValuesAndIndexArray(eigenVecs, vectorNumbers);
                    eigenVecsK = create_submatrix_Keep_cols(eigenVecs, selectedColumns);
                    eigenfaces = multiply_Matrices(phi, eigenVecsK);
                    normEigenfaces = normalize(eigenfaces);

                    weightsMatrix = new double[phi.length][phi[0].length];
                    populateWeightsMatrix(weightsMatrix, phi, normEigenfaces);

                    phiVector = centralizeVector(matrixToArray1D(matrixCSVDouble), meanVector);
                    weightsVetorPrincipal = calculateWeights(phiVector, transposed_Matrix(normEigenfaces));

                    euclideanDistances = calculate_Euclidian_Distance(weightsVetorPrincipal, weightsMatrix);
                    closestImageIndex = check_Closer_Vetor(euclideanDistances);

                    double[] reconstructedImage = reconstructImage(meanVector, normEigenfaces, getColumn(weightsMatrix, closestImageIndex));

                    System.out.println("A imagem mais próxima é: " + csvFiles[closestImageIndex]);
                    saveImage(array1DToMatrix(reconstructedImage, imageMatrices[closestImageIndex]), csvFiles[closestImageIndex], "Output/Func3");
                    break;
                default:
                    print_Header_Function("Número inválido");
                    error_General("Erro: Opção inválida");
                    break;
            }


        }

        private static void populateLinearizedImages(double[][] linearizedImages, double[][][] imageMatrices) {
            for (int img = 0; img < imageMatrices.length; img++) {
                double[] linearizedMatrix = matrixToArray1D(imageMatrices[img]);
                for (int i = 0; i < linearizedMatrix.length; i++) {
                    linearizedImages[i][img] = linearizedMatrix[i];
                }
            }
        }

        private static void populateWeightsMatrix(double[][] weightsMatrix, double[][] phi, double[][] normEigenfaces) {
            for (int img = 0; img < phi[0].length; img++) {
                double[] weights = calculateWeightsOne(getColumn(phi, img), normEigenfaces);
                for (int i = 0; i < weights.length; i++) {
                    weightsMatrix[i][img] = weights[i];
                }
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

        //* Printar Funcionalidades
        public static void print_Function_1(double[][] matrixCSVDouble, int vectorNumbers, double[][] newEigenVectorsK, double[][] newEigenValuesK, double[][] newEigenVectorsTransposeK, double[][] resultingMatrixAk, double errorAbsMed) {
            print_Matrix(matrixCSVDouble, "Matriz Original");
            print_Matrix(resultingMatrixAk, "Matriz Resultante k:" + vectorNumbers);
            print_Matrix(newEigenValuesK, "Matriz Valores Próprios k:" + vectorNumbers);
            print_Matrix(newEigenVectorsK, "Matriz Vetores Próprios k:" + vectorNumbers);
            print_Matrix(newEigenVectorsTransposeK, "Matriz Vetores Próprios Transpostos k:" + vectorNumbers);
            System.out.println("Erro Absoluto Médio: " + errorAbsMed);
        }
        public static void print_Header_Function(String functionName) {
            System.out.println();
            print_Line(1, "-------------------------------------------------------");
            System.out.println(functionName);
            print_Line(1, "-------------------------------------------------------");
            System.out.println();
        }
        //* ----------------- Fim printar funcionalidades -----------------------

        //* -------------------- Funcionalidade 1 -----------------------
        // eu (gabriel) apenas colocarei os metodos que ainda não tem no main, depois podemos organizar melhor
        public static EigenDecomposition decompose_Matrix(double[][] arrayToDecompose) {
            Array2DRowRealMatrix matrix = new Array2DRowRealMatrix(arrayToDecompose);
            return new EigenDecomposition(matrix);
        }
        public static double calculateMAE(double[][] A, double[][] Ak) {
            int M = A.length;
            int N = A[0].length;
            double errorAbsMed = 0;
            // Percorre cada elemento da matriz
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < N; j++) {
                    errorAbsMed += Math.abs(A[i][j] - Ak[i][j]);
                }
            }

            // Calcula o erro médio
            return errorAbsMed / (M * N);
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
        // TODO verificar se o metodo é mais eficiente que o método keepColumns
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
        // TODO ----------------------------------------------------------------
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

        private static void saveMatrixToFile(double[][] matrix, String inputCsvPath, String outputFolderPath) {
            File outputFolder = new File(outputFolderPath);
            if (!outputFolder.exists()) {
                if (outputFolder.mkdirs()) {
                    System.out.println("Diretório criado: " + outputFolderPath);
                } else {
                    System.err.println("Falha ao criar o diretório: " + outputFolderPath);
                    return;
                }
            }

            String csvFileName = new File(inputCsvPath).getName();
            String newFileName = "Reconstruida-" + csvFileName;
            String outputPath = outputFolderPath + "/" + newFileName;

            int counter = 1;
            File file = new File(outputPath);
            while (file.exists()) {
                file = new File(outputFolderPath + "/" + newFileName.replace(".csv", "(" + counter + ").csv"));
                counter++;
            }

            try (PrintWriter writer = new PrintWriter(file)) {
                for (double[] row : matrix) {
                    String rowString = String.join(" ; ", Arrays.stream(row)
                            .mapToObj(val -> String.format("%.3f", val))
                            .toArray(String[]::new));
                    writer.println(rowString);
                }
                System.out.println("Arquivo CSV criado com sucesso: " + file.getName());
            } catch (IOException e) {
                System.err.println("Erro ao salvar a matriz no arquivo: " + e.getMessage());
                e.printStackTrace();
            }
        }
        //* ----------------- Fim funcionalidade 1 ------------------

        //* ----------------- Funcionalidade 2 ------------------
        public static double[][][] get_Matrices_From_CSV_Folder(String folderLocation) {
            File folder = new File(folderLocation);
            File[] csvFiles = folder.listFiles((dir, name) -> name.toLowerCase().endsWith(".csv"));
            if (csvFiles == null || csvFiles.length == 0) {
                throw new RuntimeException("Nenhum arquivo CSV encontrado na pasta: " + folderLocation);
            }

            for (int i = 0; i < csvFiles.length - 1; i++) {
                for (int j = i + 1; j < csvFiles.length; j++) {
                    if (csvFiles[i].getName().compareTo(csvFiles[j].getName()) > 0) {
                        File temp = csvFiles[i];
                        csvFiles[i] = csvFiles[j];
                        csvFiles[j] = temp;
                    }
                }
            }

            double[][][] matrices = new double[csvFiles.length][][];

            for (int i = 0; i < csvFiles.length; i++) {
                File csvFile = csvFiles[i];
                matrices[i] = readCSVToMatrix(csvFile.getPath());
            }

            return matrices;
        }
        private static double[][] readCSVToArray(String filePath) {
            ArrayList<double[]> rows = new ArrayList<>();
            try (Scanner scanner = new Scanner(new File(filePath))) {
                while (scanner.hasNextLine()) {
                    String line = scanner.nextLine().trim();
                    if (!line.isEmpty()) {
                        String[] values = line.split(",");
                        double[] row = new double[values.length];
                        for (int i = 0; i < values.length; i++) {
                            row[i] = Double.parseDouble(values[i].trim());
                        }
                        rows.add(row);
                    }
                }
            } catch (FileNotFoundException e) {
                throw new RuntimeException("Erro: Arquivo CSV não encontrado no caminho: " + filePath, e);
            } catch (NumberFormatException e) {
                throw new RuntimeException("Erro: Valores inválidos encontrados no arquivo CSV.", e);
            }

            // Converte a lista de linhas em uma matriz
            double[][] matrix = new double[rows.size()][];
            for (int i = 0; i < rows.size(); i++) {
                matrix[i] = rows.get(i);
            }
            return matrix;
        }
        public static double[] centralizeVector(double[] vector, double[] meanVector) {
            if (vector.length != meanVector.length) {
                throw new IllegalArgumentException("O comprimento do vetor deve ser igual ao tamanho do vetor médio.");
            }

            double[] phi = new double[vector.length]; // Vetor para armazenar o vetor centralizado

            for (int i = 0; i < vector.length; i++) {
                phi[i] = vector[i] - meanVector[i]; // Centraliza o valor
            }

            return phi; // Retorna o vetor centralizado
        }
        public static double[] calculate_Euclidian_Distance(double[] vetorPrincipal, double[][] matrizVetores) {
            double[] resultado = new double[matrizVetores[0].length];
            for (int i = 0; i < matrizVetores[0].length; i++) {
                double soma = 0;
                for (int j = 0; j < matrizVetores.length; j++) {
                    soma += Math.pow(vetorPrincipal[j] - matrizVetores[j][i], 2);
                }
                resultado[i] = Math.sqrt(soma);
            }
            return resultado;
        }
        public static int check_Closer_Vetor(double[] resultado) {
            double min = resultado[0];
            int min_Pos = 0;
            for (int j = 1; j < resultado.length; j++) {
                if (resultado[j] < min) {
                    min = resultado[j];
                    min_Pos = j;
                }
            }
            return min_Pos;
        }
        public static double[] getColumn(double[][] matrix, int column) {
            double[] columnData = new double[matrix.length];
            for (int i = 0; i < matrix.length; i++) {
                columnData[i] = matrix[i][column];
            }
            return columnData;
        }
        public static String[] getCSVFileNames(String folderLocation) {
            File folder = new File(folderLocation);
            File[] csvFiles = folder.listFiles((dir, name) -> name.toLowerCase().endsWith(".csv"));
            if (csvFiles == null || csvFiles.length == 0) {
                throw new RuntimeException("Nenhum arquivo CSV encontrado na pasta: " + folderLocation);
            }

            String[] fileNames = new String[csvFiles.length];
            for (int i = 0; i < csvFiles.length; i++) {
                fileNames[i] = csvFiles[i].getName();
            }

            return fileNames;
        }

        public static double[] matrixToArray1D(double[][] matrix) {
            int rows = matrix.length;
            int columns = matrix[0].length;
            double[] array = new double[rows * columns];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    array[i * columns + j] = matrix[i][j];
                }
            }
            return array;
        }
        public static double[] calculateMeanVector(double[][] linearizedImages) {
            int numPixels = linearizedImages.length;
            int numImages = linearizedImages[0].length;
            double[] meanVector = new double[numPixels];

            for (int i = 0; i < numPixels; i++) {
                double sum = 0;
                for (int j = 0; j < numImages; j++) {
                    sum += linearizedImages[i][j];
                }
                meanVector[i] = sum / numImages;
            }
            return meanVector;
        }
        public static void saveImage(double[][] imageArray, String inputCsvPath, String outputFolderPath) {
            int height = imageArray.length;
            int width = imageArray[0].length;

            double min = Double.MAX_VALUE;
            double max = Double.MIN_VALUE;

            for (double[] row : imageArray) {
                for (double val : row) {
                    if (val < min) min = val;
                    if (val > max) max = val;
                }
            }

            int[][] normalizedImage = new int[height][width];
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    normalizedImage[y][x] = (int) ((imageArray[y][x] - min) / (max - min) * 255);
                }
            }

            String csvFileName = new File(inputCsvPath).getName();
            String pngFileName = csvFileName.replace(".csv", ".jpg");
            String outputPath = outputFolderPath + "/" + pngFileName;

            File outputFolder = new File(outputFolderPath);
            if (!outputFolder.exists()) {
                if (!outputFolder.mkdirs()) {
                    System.err.println("Falha ao criar o diretório: " + outputFolderPath);
                    return;
                }
            }

            int counter = 1;
            File file = new File(outputPath);
            while (file.exists()) {
                file = new File(outputFolderPath + "/" + pngFileName.replace(".png", "(" + counter + ").png"));
                counter++;
            }

            try {
                writeArrayAsImage(normalizedImage, file.getAbsolutePath());
                System.out.println("Imagem salva com sucesso: " + file.getAbsolutePath());
            } catch (IOException e) {
                System.err.println("Erro ao salvar a imagem: " + e.getMessage());
            }
        }
        public static double[][] centralizeImages(double[][] images, double[] meanVector) {
            int numPixels = meanVector.length;      // Número de pixels por imagem
            int numImages = images[0].length;      // Número de imagens

            if (images.length != numPixels) {
                throw new IllegalArgumentException("O número de linhas em 'images' deve ser igual ao tamanho do 'meanVector'.");
            }

            double[][] phi = new double[numPixels][numImages]; // Matriz para armazenar as imagens centralizadas

            for (int i = 0; i < numPixels; i++) { // Para cada pixel
                for (int j = 0; j < numImages; j++) { // Para cada imagem
                    phi[i][j] = images[i][j] - meanVector[i]; // Centraliza o valor
                }
            }

            return phi; // Retorna a matriz centralizada
        }
        public static double[] calculateWeightsOne (double[] phi, double[][] matrixU) {
            System.out.println("phi: " + phi.length);
            System.out.println("matrixU: " + matrixU.length);
            System.out.println("matrixU[0]: " + matrixU[0].length);
            if (phi.length != matrixU.length) {
                throw new IllegalArgumentException("O comprimento de 'phi' deve ser igual ao número de linhas em 'eigenfaces'.");
            }

            double[] weights = new double[matrixU[0].length];

            for (int j = 0; j < matrixU[0].length; j++) {
                for (int i = 0; i < matrixU.length; i++) {
                    weights[j] += phi[i] * matrixU[i][j];
                }
            }
            return weights;
        }
        public static double[] calculateWeights(double[] phi, double[][] matrixU) {
            if (phi.length != matrixU[0].length) {
                throw new IllegalArgumentException("O comprimento de 'phi' deve ser igual ao número de linhas em 'eigenfaces'.");
            }

            double[] weights = new double[matrixU[0].length];

            for (int i = 0; i < matrixU.length; i++) {
                for (int j = 0; j < matrixU[0].length; j++) {
                    weights[i] += phi[j] * matrixU[i][j];
                }
            }
            return weights;
        }
        public static double[] reconstructImage(double[] meanVector, double[][] eigenfaces, double[] weights) {
            int quantityEigenfaces = eigenfaces[0].length;

            double[] reconstructed = new double[meanVector.length];
            for (int i = 0; i < meanVector.length; i++) {
                reconstructed[i] = meanVector[i];
            }

            for (int j = 0; j < quantityEigenfaces; j++) {
                for (int i = 0; i < eigenfaces.length; i++) {
                    reconstructed[i] += weights[j] * eigenfaces[i][j];
                }
            }

            return reconstructed;
        }
        public static double[][] array1DToMatrix(double[] array, double[][] originalMatrix) {
            int rows = originalMatrix.length;
            int columns = originalMatrix[0].length;
            if (array.length != rows * columns) {
                throw new IllegalArgumentException("O tamanho do vetor não corresponde às dimensões da matriz.");
            }
            double[][] matrix = new double[rows][columns];
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    matrix[i][j] = array[i * columns + j];
                }
            }
            return matrix;
        }
        public static void writeArrayAsImage(int[][] array, String outputFilePath) throws IOException {
            int height = array.length;
            int width = array[0].length;

            // Create a BufferedImage
            BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);

            // Set the pixel intensities
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    int intensity = array[y][x];
                    if (intensity < 0 || intensity > 255) {
                        throw new IllegalArgumentException("Pixel intensity must be between 0 and 255.");
                    }
                    int rgb = (intensity << 16) | (intensity << 8) | intensity; // Set the same value for R, G, B
                    image.setRGB(x, y, rgb);
                }
            }

            // Write the image to the file
            File outputFile = new File(outputFilePath);
            ImageIO.write(image, "png", outputFile);
        }
        private static double[][] readCSVToMatrix(String path) {
            try {
                // Conta o número de linhas no arquivo
                Scanner lineCounter = new Scanner(new File(path));
                int rowCount = 0;
                int columnCount = 0;

                while (lineCounter.hasNextLine()) {
                    String line = lineCounter.nextLine();
                    if (!line.trim().isEmpty()) {
                        rowCount++;
                        if (columnCount == 0) {
                            columnCount = line.split(",").length;
                        }
                    }
                }
                lineCounter.close();

                // Inicializa a matriz
                double[][] matrix = new double[rowCount][columnCount];

                // Lê o arquivo novamente para preencher a matriz
                Scanner fileScanner = new Scanner(new File(path));
                int row = 0;
                while (fileScanner.hasNextLine()) {
                    String line = fileScanner.nextLine();
                    if (!line.trim().isEmpty()) {
                        String[] values = line.split(",");
                        for (int col = 0; col < values.length; col++) {
                            matrix[row][col] = Double.parseDouble(values[col].trim());
                        }
                        row++;
                    }
                }
                fileScanner.close();
                return matrix;
            } catch (Exception e) {
                throw new RuntimeException("Erro ao ler o arquivo CSV: " + e.getMessage(), e);
            }
        }
        public static double[][] covariances(double[][] matrixA) {
            int quantityOfImages = matrixA[0].length;
            double[][] matrixATransposed = transposed_Matrix(matrixA);;
            double[][] matrixATmultiplyByA = multiply_Matrices(matrixATransposed, matrixA);
            return multiply_Matrix_Escalar(matrixATmultiplyByA, 1.0 / quantityOfImages);
        }
        public static double[][] eigenVectors(double[][] matrix) {
            double[][] matrixTransposed = transposed_Matrix(matrix);
            double[][] eigenVectors = multiply_Matrices(matrixTransposed, matrix);

            //! Verificar se o decompose funciona
            EigenDecomposition decomposedMatrixThree = decompose_Matrix(eigenVectors);
            RealMatrix eigenVectorsMatrix = decomposedMatrixThree.getV();
            return eigenVectorsMatrix.getData();
        }
        public static double[][] normalize(double[][] eigenVectorsATxA) {
            for (int i = 0; i < eigenVectorsATxA[0].length; i++) {
                double norm = 0;

                for (int j = 0; j < eigenVectorsATxA.length; j++) {
                    norm += eigenVectorsATxA[j][i] * eigenVectorsATxA[j][i];
                }
                norm = Math.sqrt(norm);

                for (int j = 0; j < eigenVectorsATxA.length; j++) {
                    eigenVectorsATxA[j][i] /= norm;
                }
            }
            return eigenVectorsATxA;
        }

        //* ----------------- Fim funcionalidade 2 ------------------


        //! ------------------ Error Messages ------------------
        public static void error_General(String error) {
            System.out.println(error);
            System.exit(1);
        }
        //! ------------------ Fim error messages --------------
    }