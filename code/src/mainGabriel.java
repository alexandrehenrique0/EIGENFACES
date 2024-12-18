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

public class mainGabriel {

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
        int vectorNumbers = -1;
        String csvLocation = null;
        String imageFolderLocation = null;

        // VERIFICAÇÃO DE PARÂMETROS (ARGS) PARA MÉTODOS DE INPUT
        if (checkCorrectParametersStructure(args)) {

            // Receber os parâmetros
            function = receiveFunction(args);
            vectorNumbers = receiveNumberVectors(args);
            csvLocation = receiveCsvLocation(args);
            imageFolderLocation = receiveImageLocation(args);

            // Verificar se os arquivo e diretório existem
            checkExistanceFileDirectory(csvLocation, imageFolderLocation);

            // Obter a matriz do CSV
            double[][] matrixCSVDouble = getMatrixFromCsv(csvLocation);

            // Obter a matriz do CSV para a função 2
            double[][][] matrixCSVDouble3D = getMatricesFromCsvFolder(imageFolderLocation);

            // Função que contém as funções principais
            switchPrimaryFunctions(function, matrixCSVDouble, vectorNumbers, csvLocation, matrixCSVDouble3D, imageFolderLocation);

        } else if (args.length == 0) {
            // Mostrar as opções num menu e receber os parâmetros
            do {
                uiFunctionParameterMenu();
                function = receiveFunction(null);
            } while (!check_function(function));

            if (function == 4) {
                uiQuitParameterMenu();
                receiveExitConfirmation(null);
            }

            uiVectorNumbersParameterMenu();
            vectorNumbers = receiveNumberVectors(null);

            uiCsvLocationParameterMenu();
            csvLocation = receiveCsvLocation(null);

            uiImageLocationParameterMenu();
            imageFolderLocation = receiveImageLocation(null);

            // Verificar se os arquivo e diretório existem
            checkExistanceFileDirectory(csvLocation, imageFolderLocation);

            // Obter a matriz do CSV
            double[][] matrixCSVDouble = getMatrixFromCsv(csvLocation);

            // Obter a matriz do CSV para a função 2
            double[][][] matrixCSVDouble3D = getMatricesFromCsvFolder(imageFolderLocation);

            // Função que contém as funções principais
            switchPrimaryFunctions(function, matrixCSVDouble, vectorNumbers, csvLocation, matrixCSVDouble3D, imageFolderLocation);


        } else {
            errorGeneral("Erro: Parâmetros inválidos");
        }
    }

    public static void switchPrimaryFunctions(int function, double[][] matrixCSVDouble, int vectorNumbers, String csvLocation, double[][][] matrixCSVDouble3D, String imageFolderLocation) {
        // Common variables for all 3 functions
        String[] csvFiles = getCSVFileNames(imageFolderLocation);
        double[][] linearizedImages = new double[matrixCSVDouble3D[0].length * matrixCSVDouble3D[0].length][matrixCSVDouble3D.length];
        double[] averageVectors;
        double[][] phi;
        double[][] covariance;
        double[][] eigenVectors;
        double[][] selectedColumnsK;
        double[][] newEigenVectorsK;
        double[][] expandedVectorsK;
        double[][] eigenfaces;
        double[][] weightsMatrix = new double[matrixCSVDouble3D.length][matrixCSVDouble3D[0].length * matrixCSVDouble3D[0].length];
        double[] phiVector;
        double[] principalWeightsVector;
        double[] euclideanDistances;
        int closestImageIndex;

        populateLinearizedImages(linearizedImages, matrixCSVDouble3D);

        switch (function) {
            case 1:
                printHeaderFunction("Decomposição Própria de uma Matriz Simétrica");
                eigenVectors = getEigenVectors(matrixCSVDouble);
                double[][] eigenValues = getEigenValues(matrixCSVDouble);
                double[][] valuesAndIndexArray = getValuesAndIndexArray(eigenValues, vectorNumbers);

                newEigenVectorsK = createSubmatrix(eigenVectors, valuesAndIndexArray);
                double[][] newEigenValuesK = constructDiagonalMatrix(valuesAndIndexArray);
                double[][] newEigenVectorsTransposeK = transposeMatrix(newEigenVectorsK);
                double[][] matrixEigenFaces = multiplyVectorsValuesVectorsTransposed(newEigenVectorsK, newEigenValuesK, newEigenVectorsTransposeK);

                double maximumAbsolutError = calculateMAE(matrixCSVDouble, matrixEigenFaces);

                saveMatrixToFile(matrixEigenFaces, csvLocation, "Output/Func1");
                printFunction1(matrixCSVDouble, vectorNumbers, newEigenVectorsK, newEigenValuesK, newEigenVectorsTransposeK, matrixEigenFaces, maximumAbsolutError);
                break;
            case 2:
                printHeaderFunction("Reconstrução de Imagens usando Eigenfaces");

                averageVectors = calculateAverageVector(linearizedImages);
                phi = centralizeImages(linearizedImages, averageVectors);
                double[][] phiT = transposeMatrix(phi);
                double[][] phiTxPhi = multiplyMatrices(phiT, phi);
                eigenVectors = getEigenVectors(phiTxPhi);

                if (vectorNumbers == -1 || vectorNumbers > eigenVectors[0].length) {
                    vectorNumbers = eigenVectors[0].length;
                }

                selectedColumnsK = getValuesAndIndexArray(eigenVectors, vectorNumbers);
                newEigenVectorsK = createSubmatrix(eigenVectors, selectedColumnsK);
                expandedVectorsK = multiplyMatrices(phi, newEigenVectorsK);
                eigenfaces = normalize(expandedVectorsK);

                populateWeightsMatrix(weightsMatrix, phi, eigenfaces);

                for (int img = 0; img < linearizedImages[0].length; img++) {
                    double[] columnWeights = getColumn(weightsMatrix, img);
                    double[] reconstructedImage = reconstructImage(averageVectors, eigenfaces, columnWeights);
                    double[][] reconstructedImageMatrix = array1DToMatrix(reconstructedImage, matrixCSVDouble3D[img]);

                    saveImage(reconstructedImageMatrix, csvFiles[img], "Output/Func2/ImagensReconstruidas");
                    saveMatrixToFile(reconstructedImageMatrix, csvFiles[img], "Output/Func2/Eigenfaces");
                }
                break;
            case 3:
                printHeaderFunction("Identificação de imagem mais próxima");

                averageVectors = calculateAverageVector(linearizedImages);
                phi = centralizeImages(linearizedImages, averageVectors);
                covariance = covariances(phi);
                eigenVectors = getEigenVectors(covariance);

                if (vectorNumbers == -1 || vectorNumbers > eigenVectors[0].length) {
                    vectorNumbers = eigenVectors[0].length;
                }

                selectedColumnsK = getValuesAndIndexArray(eigenVectors, vectorNumbers);
                newEigenVectorsK = createSubmatrix(eigenVectors, selectedColumnsK);
                expandedVectorsK = multiplyMatrices(phi, newEigenVectorsK);
                eigenfaces = normalize(expandedVectorsK);
                populateWeightsMatrix(weightsMatrix, phi, eigenfaces);

                double[] matrixLinearized = matrixToArray1D(matrixCSVDouble);
                phiVector = centralizeVector(matrixLinearized, averageVectors);

                double[][] eigenfacesTransposed = transposeMatrix(eigenfaces);
                principalWeightsVector = calculateWeights(phiVector, eigenfacesTransposed);

                euclideanDistances = calculateEuclidianDistance(principalWeightsVector, weightsMatrix);
                closestImageIndex = checkCloserVetor(euclideanDistances);

                double[] columnWeights = getColumn(weightsMatrix, closestImageIndex);

                double[] reconstructedImage = reconstructImage(averageVectors, eigenfaces, columnWeights);

                System.out.println("A imagem mais próxima é: " + csvFiles[closestImageIndex]);

                double[][] reconstructedImageMatrix = array1DToMatrix(reconstructedImage, matrixCSVDouble3D[0]);

                saveImage(reconstructedImageMatrix, csvFiles[closestImageIndex], "Output/Func3");
                break;
            case 4:
                printHeaderFunction("Tem certeza ?");
                break;
            default:
                printHeaderFunction("Número inválido");
                System.out.println("Por favor, insira uma opção válida.");
                uiFunctionParameterMenu();
                receiveFunction(null);
                switchPrimaryFunctions(function, matrixCSVDouble, vectorNumbers, csvLocation, matrixCSVDouble3D, imageFolderLocation);
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

    private static void populateWeightsMatrix(double[][] weightsMatrix, double[][] phi, double[][] eigenfaces) {
        for (int img = 0; img < phi[0].length; img++) {
            double[] weights = calculateWeightsOne(getColumn(phi, img), eigenfaces);
            for (int i = 0; i < weights.length; i++) {
                weightsMatrix[i][img] = weights[i];
            }
        }
    }

    //* ------------------ Verificações ------------------
    public static boolean checkCorrectParametersStructure(String[] parameters) {
        if (parameters.length == 8) {
            return parameters[0].equals("-f") && parameters[2].equals("-k") && parameters[4].equals("-i") && parameters[6].equals("-j");
        }
        return false;
    }

    public static boolean check_function(int function) {
        return function >= 1 && function <= 4;
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

    public static void checkExistanceFileDirectory(String csvLocation, String imageFolderLocation) {
        try {
            SCANNER_CSV = new Scanner(new File(csvLocation));
            File[] IMAGE_FILES = new File(imageFolderLocation).listFiles();
        } catch (FileNotFoundException e) {
            errorGeneral("Erro ao abrir os arquivos: " + e.getMessage());
        }
    }
    //* ------------------ Fim verificações ------------------


    //* ------------------ Menus de opções ------------------
    public static void uiFunctionParameterMenu() {
        System.out.println("------------- Que função deseja realizar? -------------");
        System.out.println("1 - Decomposição Própria de uma Matriz Simétrica");
        System.out.println("2 - Reconstrução de Imagens usando Eigenfaces");
        System.out.println("3 - Identificação de imagem mais próxima");
        System.out.println("4 - Deseja sair da aplicação ? ");
        System.out.println("-------------------------------------------------------");
        System.out.printf("Opção: ");
    }

    public static void uiVectorNumbersParameterMenu() {
        System.out.println("------ Quantos vetores próprios deseja utilizar? ------");
        System.out.printf("Quantidade: ");
    }

    public static void uiCsvLocationParameterMenu() {
        System.out.println("---- Qual a localização do csv que deseja utilizar? ---");
        System.out.printf("Localização: ");
    }

    public static void uiImageLocationParameterMenu() {
        System.out.println("-------- Qual a localização da base de imagens? -------");
        System.out.printf("Localização: ");
    }

    public static void uiQuitParameterMenu() {
        System.out.println("-------------------------------------------------------");
        System.out.println("Deseja sair da aplicação? (S/N)");
        System.out.printf("Opção: ");
    }

    //* --------------------- Fim menus de opções ------------------


    //* ------------------ Receber parâmetros ------------------
    public static int receiveFunction(String[] args) {
        int functionArgs;
        if (args == null) {
            int function = SCANNER.nextInt();
            if (!check_function(function)) {
                System.out.println("Erro: Opção inválida.");
                System.out.println("Retornando ao menu de opções.");
                System.out.println();
            }
            return function;
        } else {
            functionArgs = Integer.parseInt(args[1]);
            if (!check_function(functionArgs)) {
                errorGeneral("Erro: Opção inválida");
            }
            return functionArgs;
        }
    }

    public static int receiveNumberVectors(String[] args) {
        int vectorNumbersArgs;
        if (args == null) {
            return SCANNER.nextInt();
        } else {
            vectorNumbersArgs = Integer.parseInt(args[3]);
            return vectorNumbersArgs;
        }
    }

    public static String receiveCsvLocation(String[] args) {
        String csvLocationArgs;
        if (args == null) {
            String csvLocation = SCANNER.next();
            if (!check_csvLocation(csvLocation)) {
                errorGeneral("Erro: Localização inválida csv");
            }
            return csvLocation;
        } else {
            csvLocationArgs = args[5];
            if (!check_csvLocation(csvLocationArgs)) {
                errorGeneral("Erro: Localização inválida csv");
            }
            return csvLocationArgs;
        }
    }

    public static String receiveImageLocation(String[] args) {
        String imageFolderLocationArgs;
        if (args == null) {
            String imageFolderLocation = SCANNER.next();
            if (!check_imageFolderLocation(imageFolderLocation)) {
                errorGeneral("Erro: Localização inválida da pasta de imagens");
            }
            return imageFolderLocation;
        } else {
            imageFolderLocationArgs = args[7];
            if (!check_imageFolderLocation(imageFolderLocationArgs)) {
                errorGeneral("Erro: Localização inválida da pasta de imagens");
            }
            return imageFolderLocationArgs;
        }
    }

    public static void receiveExitConfirmation(String[] args) {
        if (args == null) {
            String confirmeExit;
            do {
                confirmeExit = SCANNER.next().toUpperCase();
                if (confirmeExit.equals("S")) {
                    System.exit(0);
                } else if (confirmeExit.equals("N")) {
                    System.out.println("Retornando ao menu de opções.");
                    uiFunctionParameterMenu();
                } else {
                    System.out.println("Erro: Responda com S/N");
                }
            } while (!confirmeExit.equals("N"));
        }
    }
    //* -------------------- Fim receber parâmetros ------------------


    //* ------------------ Leitura de CSV ------------------
    public static double[][] getMatrixFromCsv(String csvLocation) {
        int[] dimensions = get_Dimensions();
        int rows = dimensions[0];
        int cols = dimensions[1];

        double[][] matrix = new double[rows][cols];
        populateMatrix(matrix, csvLocation);

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
            errorGeneral("Erro: Dimensões da matriz fora dos limites: " + rows + "x" + cols);
        }
        return new int[]{rows, cols};
    }

    private static boolean check_Size_Boundaries(int rows, int cols) {
        return rows > MAX_SIZE_ROWS || cols > MAX_SIZE_COLS || rows < MIN_SIZE_ROWS || cols < MIN_SIZE_COLS;
    }

    private static void populateMatrix(double[][] matrix, String csvLocation) {
        try {
            SCANNER_CSV = new Scanner(new File(csvLocation));
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Erro ao reabrir o arquivo CSV: " + e.getMessage());
        }
        int row = 0;
        while (SCANNER_CSV.hasNextLine()) {
            String line = SCANNER_CSV.nextLine().trim();
            if (!line.isEmpty()) {
                populateRow(matrix, row, line);
                row++;
            }
        }
        SCANNER_CSV.close();
    }

    private static void populateRow(double[][] matrix, int row, String line) {
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
    public static double[][] sumMatrices(double[][] matrizLeft, double[][] matrizRight) {
        double[][] matrizResultante = new double[matrizLeft.length][matrizLeft[0].length];
        for (int i = 0; i < matrizLeft.length; i++) {
            for (int j = 0; j < matrizLeft[0].length; j++) {
                matrizResultante[i][j] = matrizLeft[i][j] + matrizRight[i][j];
            }
        }
        return matrizResultante;
    }

    public static double[][] subtractionMatrices(double[][] matrizLeft, double[][] matrizRight) {
        double[][] matrizResultante = new double[matrizLeft.length][matrizLeft[0].length];
        for (int i = 0; i < matrizLeft.length; i++) {
            for (int j = 0; j < matrizLeft[0].length; j++) {
                matrizResultante[i][j] = matrizLeft[i][j] - matrizRight[i][j];
            }
        }
        return matrizResultante;
    }

    public static double[][] multiplyMatrices(double[][] matrizLeft, double[][] matrizRight) {
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

    public static double[][] multiplyMatrixEscalar(double[][] matriz, double escalar) {
        double[][] matrizResultante = new double[matriz.length][matriz[0].length];
        for (int i = 0; i < matriz.length; i++) {
            for (int j = 0; j < matriz[0].length; j++) {
                matrizResultante[i][j] = matriz[i][j] * escalar;
            }
        }
        return matrizResultante;
    }

    public static double[][] transposeMatrix(double[][] matrix) {
        double[][] transposedMatrix = new double[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                transposedMatrix[j][i] = matrix[i][j];
            }
        }
        return transposedMatrix;
    }

    public static double[][] createSubmatrix(double[][] eigenVectors, double[][] valuesAndIndexArray) {

        boolean[] keepColumnsBOL = new boolean[eigenVectors[0].length];

        for (double[] col : valuesAndIndexArray) {
            keepColumnsBOL[(int) col[1]] = true;
        }

        double[][] submatrix = new double[eigenVectors.length][valuesAndIndexArray.length];

        int sub_i = 0;
        for (double[] doubles : eigenVectors) {
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
    public static void printMatrix(double[][] matrixToPrint, String matrixName) {
        System.out.println("Matriz: " + matrixName + " ↓");
        printLine(matrixToPrint[0].length, "____________");

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
        printLine(matrixToPrint[0].length, "============");
        System.out.println();
    }

    // Metodo para printar linhas de caracteres no console.
    private static void printLine(int length, String pattern) {
        for (int i = 0; i < length; i++) {
            System.out.print(pattern);
        }
        System.out.println();
    }
    //* ----------------- Fim printar matrizes -----------------------

    //* Printar Funcionalidades
    public static void printFunction1(double[][] matrixCSVDouble, int vectorNumbers, double[][] newEigenVectorsK, double[][] newEigenValuesK, double[][] newEigenVectorsTransposeK, double[][] resultingMatrixAk, double errorAbsMed) {
        printMatrix(matrixCSVDouble, "Matriz Original");
        printMatrix(resultingMatrixAk, "Matriz Resultante k:" + vectorNumbers);
        printMatrix(newEigenValuesK, "Matriz Valores Próprios k:" + vectorNumbers);
        printMatrix(newEigenVectorsK, "Matriz Vetores Próprios k:" + vectorNumbers);
        printMatrix(newEigenVectorsTransposeK, "Matriz Vetores Próprios Transpostos k:" + vectorNumbers);
        System.out.println("Erro Absoluto Médio: " + errorAbsMed);
    }

    public static void printHeaderFunction(String functionName) {
        System.out.println();
        printLine(1, "-------------------------------------------------------");
        System.out.println(functionName);
        printLine(1, "-------------------------------------------------------");
        System.out.println();
    }
    //* ----------------- Fim printar funcionalidades -----------------------

    //* -------------------- Funcionalidade 1 -----------------------
    // eu (gabriel) apenas colocarei os metodos que ainda não tem no main, depois podemos organizar melhor
    public static EigenDecomposition decomposeMatrix(double[][] arrayToDecompose) {
        Array2DRowRealMatrix matrix = new Array2DRowRealMatrix(arrayToDecompose);
        return new EigenDecomposition(matrix);
    }

    public static double calculateMAE(double[][] originalMatrix, double[][] matrixEigenFaces) {
        int M = originalMatrix.length;
        int N = originalMatrix[0].length;
        double errorAbsMed = 0;
        // Percorre cada elemento da matriz
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                errorAbsMed += Math.abs(originalMatrix[i][j] - matrixEigenFaces[i][j]);
            }
        }

        // Calcula o erro médio
        return errorAbsMed / (M * N);
    }

    public static double[][] getEigenVectors(double[][] matrix) {
        EigenDecomposition eigenDecomposition = decomposeMatrix(matrix);
        RealMatrix eigenVectors = eigenDecomposition.getV();
        return eigenVectors.getData();
    }

    public static double[][] getEigenValues(double[][] matrix) {
        EigenDecomposition eigenDecomposition = decomposeMatrix(matrix);
        RealMatrix eigenValues = eigenDecomposition.getD();
        return eigenValues.getData();
    }

    public static double[][] getEigenVectorsTranspose(double[][] matrix) {
        EigenDecomposition eigenDecomposition = decomposeMatrix(matrix);
        RealMatrix eigenVectorsTranspose = eigenDecomposition.getVT();
        return eigenVectorsTranspose.getData();
    }

    private static double[][] constructDiagonalMatrix(double[][] matrixvaluesK) {
        double[][] matrixvaluesKPrint = new double[matrixvaluesK.length][matrixvaluesK.length];
        for (int i = 0; i < matrixvaluesK.length; i++) {
            matrixvaluesKPrint[i][i] = matrixvaluesK[i][0];
        }
        return matrixvaluesKPrint;
    }

    private static double[][] getValuesAndIndexArray(double[][] eigenValuesArray, int eigenfaces) {
        double[][] valuesAndIndexArray = new double[eigenfaces][2];

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

    public static double[][] multiplyVectorsValuesVectorsTransposed(double[][] matrixVectors, double[][] matrixValues, double[][] matrixVectorsTranspose) {
        return multiplyMatrices(multiplyMatrices(matrixVectors, matrixValues), matrixVectorsTranspose);
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
    public static double[][][] getMatricesFromCsvFolder(String folderLocation) {
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

    public static double[] calculateEuclidianDistance(double[] vetorPrincipal, double[][] matrizVetores) {
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

    public static int checkCloserVetor(double[] resultado) {
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

    public static double[] calculateAverageVector(double[][] linearizedImages) {
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

    public static double[] calculateWeightsOne(double[] phi, double[][] matrixU) {
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

    public static double[] reconstructImage(double[] averageVector, double[][] eigenfaces, double[] columnWeights) {
        int quantityEigenfaces = eigenfaces[0].length;

        double[] reconstructed = new double[averageVector.length];
        for (int i = 0; i < averageVector.length; i++) {
            reconstructed[i] = averageVector[i];
        }

        for (int j = 0; j < quantityEigenfaces; j++) {
            for (int i = 0; i < eigenfaces.length; i++) {
                reconstructed[i] += columnWeights[j] * eigenfaces[i][j];
            }
        }
        return reconstructed;
    }

    public static double[][] array1DToMatrix(double[] reconstructedImage, double[][] matrixSizeExample) {
        int rows = matrixSizeExample.length;
        int columns = matrixSizeExample[0].length;
        if (reconstructedImage.length != rows * columns) {
            throw new IllegalArgumentException("O tamanho do vetor não corresponde às dimensões da matriz.");
        }
        double[][] reconstructedMatrix = new double[rows][columns];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                reconstructedMatrix[i][j] = reconstructedImage[i * columns + j];
            }
        }
        return reconstructedMatrix;
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
        double[][] matrixATransposed = transposeMatrix(matrixA);
        ;
        double[][] matrixATmultiplyByA = multiplyMatrices(matrixATransposed, matrixA);
        return multiplyMatrixEscalar(matrixATmultiplyByA, 1.0 / quantityOfImages);
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
    public static void errorGeneral(String error) {
        System.out.println(error);
        System.exit(1);
    }


    //! ------------------ Fim error messages --------------
}