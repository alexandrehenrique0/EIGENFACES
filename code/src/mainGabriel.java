import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.*;
//* NOTE: The code has been divided into several parts to facilitate reading and understanding.
//* Collapse all for a better view.
//* Use the keyboard shortcut Ctrl + F and search for "//*" to navigate between parts of the code.


public class mainGabriel {
    //* Constantes para limites de tamanho.
    public static final int MAX_SIZE_ROWS_AND_COLS = 256;
    public static final int MIN_SIZE_ROWS_AND_COLS = 1;
    public static final int MIN_QUANTITY_VECTORS = 1;
    private static final int MIN_BIT_VALUE = 0;
    private static final int MAX_BIT_VALUE = 255;

    //* Scanner global para ser utilizado em todos os métodos necessários.
    public static Scanner scanner = new Scanner(System.in);

    public static void main(String[] args) {
        // Verificação de parâmetros para decidir entre interativo e não interativo.
        if (checkCorrectParametersStructure(args)) {
            runNonInteractive(args);
        } else {
            runInterative();
        }
        scanner.close();
    }

    //* ------------------ Modos de execução ------------------
    public static void runInterative() {
        // Parâmetro de entrada
        int function;

        // Roda enquanto a função for inválida
        function = verifyFunction();

        switch (function) {
            case 1:
                function1();
                break;
            case 2:
                function2();
                break;
            case 3:
                function3();
                break;
            case 4:
                function4();
                break;
            case 5:
                devTeam();
                break;
            case 0:
                quitApplication();
                break;
        }
    }

    public static void runNonInteractive(String[] args) {
        // Parametros de entrada
        int function;
        int vectorNumbers;
        String csvLocation;
        String dataBaseLocation;

        // Receber os parâmetros
        function = receiveFunction(args);
        vectorNumbers = receiveNumberVectors(args);
        csvLocation = receiveCsvLocation(args);
        dataBaseLocation = receiveDataBaseLocation(args);

        // Verificar se os arquivo e diretório existem
        checkExistanceFileDirectory(csvLocation);

        String[] csvFiles = getCSVFileNames(dataBaseLocation);

        // Obter a matriz do CSV
        double[][] oneMatrixCsv = readCSVToMatrix(csvLocation);

        // Obter a matriz do CSV para a função 2
        double[][][] allMatricesCsv = getMatricesFromCsvFolder(dataBaseLocation);

        try {
            // Defina o caminho do arquivo com base na função escolhida
            String filePath = "Output/NaoInterativo/Func" + function;

            // Crie ou obtenha o arquivo para onde redirecionar a saída
            File file = new File(filePath, "/outputFunc" + function + ".txt");

            // Crie um PrintStream para o arquivo
            PrintStream fileOut = new PrintStream(file);

            // Redirecione System.out para o arquivo
            System.setOut(fileOut);

            // Opcional: redirecione System.err também, se necessário
            System.setErr(fileOut);

        } catch (FileNotFoundException e) {
            e.printStackTrace(); // Se o arquivo não puder ser criado, imprime erro
        }

        // Função que contém as funções principais
        runNonInterativeOutputs(function, vectorNumbers, csvLocation, csvFiles, oneMatrixCsv, allMatricesCsv);
    }
    //* ------------------ Fim modos de execução ------------------


    //* ------------------ Métodos principais ------------------
    public static void runNonInterativeOutputs(int function, int vectorNumbers, String csvLocation, String[] csvFiles, double[][] oneMatrixCsv, double[][][] allMatricesCsv) {
        double[][] linearizedImages = new double[allMatricesCsv[0].length * allMatricesCsv[0].length][allMatricesCsv.length];
        double[][] weightsMatrix = new double[allMatricesCsv[0].length * allMatricesCsv[0].length][allMatricesCsv.length];
        populateLinearizedImages(linearizedImages, allMatricesCsv);
        double[] averageVectors = calculateMeanVector(linearizedImages);
        double[][] phi = centralizeImages(linearizedImages, averageVectors);
        int vectorK = validateEigenVectors(linearizedImages, vectorNumbers);

        double[][] phiT = transposeMatrix(phi);
        double[][] phiTxPhi = multiplyMatrices(phiT, phi);
        double[][] eigenVectors = getEigenVectors(phiTxPhi);
        double[][] selectedColumnsK = getValuesAndIndexArray(eigenVectors, vectorK);
        double[][] newEigenVectorsK = createSubmatrix(eigenVectors, selectedColumnsK);
        double[][] expandedVectorsK = multiplyMatrices(phi, newEigenVectorsK);
        double[][] eigenfaces = normalize(expandedVectorsK);

        populateWeightsMatrix(weightsMatrix, phi, eigenfaces);

        switch (function) {
            case 1:
                printHeaderFunction("Decomposição Própria de uma Matriz Simétrica");
                decomposeSymmetricMatrix(vectorNumbers, csvLocation);
                System.out.println();
                System.out.println("Funcionalidade 1 finalizada.");
                break;
            case 2:
                printHeaderFunction("Reconstrução de Imagens usando Eigenfaces");
                reconstructImagesWithEigenfaces(vectorK, csvFiles, averageVectors, eigenfaces, linearizedImages, weightsMatrix, allMatricesCsv);
                System.out.println();
                System.out.println("Funcionalidade 2 finalizada.");
                break;
            case 3:
                printHeaderFunction("Identificação de imagem mais próxima");
                identifyClosestImage(vectorK, csvFiles, averageVectors, eigenfaces, oneMatrixCsv, weightsMatrix, allMatricesCsv);
                System.out.println();
                System.out.println("Funcionalidade 3 finalizada.");
                break;
        }
    }

    public static void function1() {
        int vectorNumbers = verifyVectorNumbers();
        String csvLocation = verifyCsvLocation();

        printHeaderFunction("Decomposição Própria de uma Matriz Simétrica:");

        decomposeSymmetricMatrix(vectorNumbers, csvLocation);

        System.out.println();
        System.out.println("Funcionalidade 1 finalizada, a retornar ao menu inicial.");

        runInterative();
    }

    public static void function2() {
        int vectorNumbers = verifyVectorNumbers();
        String dataBase = verifyDataBaseLocation();

        printHeaderFunction("Reconstrução de Imagens utilizando Eigenfaces");

        calculateFunction2(vectorNumbers, dataBase);

        System.out.println();
        System.out.println("Funcionalidade 2 finalizada, a retornar ao menu inicial.");

        runInterative();
    }

    public static void function3() {
        int vectorNumbers = verifyVectorNumbers();
        String csvLocation = verifyCsvLocation();
        String dataBase = verifyDataBaseLocation();

        printHeaderFunction("Identificação da imagem mais próxima utilizando Eigenfaces");

        calculateFunction3(vectorNumbers, csvLocation, dataBase);

        System.out.println();
        System.out.println("Funcionalidade 3 finalizada, a retornar ao menu inicial.");

        runInterative();
    }

    public static void function4() {
        int vectorNumbers = verifyVectorNumbers();
        String dataBase = verifyDataBaseLocation();

        printHeaderFunction("Gerar uma imagem aleatória com Eigenfaces");

        generateNewImage(vectorNumbers, dataBase);

        System.out.println();
        System.out.println("Funcionalidade 4 finalizada, a retornar ao menu inicial.");

        runInterative();
    }

    public static void devTeam(){
        printHeaderFunction("Desenvolvido por: TechTitans!");

        uiDevTeam();

        System.out.println();
        System.out.println("A retornar ao menu inicial.");

        runInterative();
    }
    //* ------------------ fim dos métodos principais ------------------


    //* ------------------ Métodos de distribuição de tarefas ------------------
    public static void decomposeSymmetricMatrix(int vectorNumbers, String csvLocation) {
        double[][] oneMatrixCsv = readCSVToMatrix(csvLocation);

        double[][] eigenVectors = getEigenVectors(oneMatrixCsv);
        double[][] eigenValues = getEigenValues(oneMatrixCsv);

        int vectorK = validateEigenVectors(oneMatrixCsv, vectorNumbers);

        double[][] valuesAndIndexArray = getValuesAndIndexArray(eigenValues, vectorK);
        double[][] newEigenVectorsK = createSubmatrix(eigenVectors, valuesAndIndexArray);
        double[][] newEigenValuesK = constructDiagonalMatrix(valuesAndIndexArray);
        double[][] newEigenVectorsTransposeK = transposeMatrix(newEigenVectorsK);
        double[][] matrixEigenFaces = multiplyMatrices(multiplyMatrices(newEigenVectorsK, newEigenValuesK), newEigenVectorsTransposeK);

        double maximumAbsolutError = calculateMAE(oneMatrixCsv, matrixEigenFaces);

        saveMatrixToFile(matrixEigenFaces, csvLocation, "Output/Func1", 1);
        printFunction1(vectorK, newEigenValuesK, newEigenVectorsK, maximumAbsolutError);
    }

    public static void calculateFunction2(int vectorNumbers, String dataBase) {
        String[] csvFiles = getCSVFileNames(dataBase);
        double[][][] allMatricesCsv = getMatricesFromCsvFolder(dataBase);

        double[][] linearizedImages = new double[allMatricesCsv[0].length * allMatricesCsv[0].length][allMatricesCsv.length];
        double[][] weightsMatrix = new double[allMatricesCsv[0].length * allMatricesCsv[0].length][allMatricesCsv.length];
        populateLinearizedImages(linearizedImages, allMatricesCsv);
        double[] averageVectors = calculateMeanVector(linearizedImages);
        double[][] phi = centralizeImages(linearizedImages, averageVectors);
        int vectorK = validateEigenVectors(linearizedImages, vectorNumbers);

        double[][] phiT = transposeMatrix(phi);
        double[][] phiTxPhi = multiplyMatrices(phiT, phi);
        double[][] eigenVectors = getEigenVectors(phiTxPhi);
        double[][] selectedColumnsK = getValuesAndIndexArray(eigenVectors, vectorK);
        double[][] newEigenVectorsK = createSubmatrix(eigenVectors, selectedColumnsK);
        double[][] expandedVectorsK = multiplyMatrices(phi, newEigenVectorsK);
        double[][] eigenfaces = normalize(expandedVectorsK);

        populateWeightsMatrix(weightsMatrix, phi, eigenfaces);

        reconstructImagesWithEigenfaces(vectorK, csvFiles, averageVectors, eigenfaces, linearizedImages, weightsMatrix, allMatricesCsv);
    }

    public static void reconstructImagesWithEigenfaces(int vectorNumbers, String[] csvFiles, double[] averageVectors, double[][] eigenfaces, double[][] linearizedImages, double[][] weightsMatrix, double[][][] allMatricesCsv) {

        System.out.println("Valores do vetor médio: " + Arrays.toString(averageVectors));
        System.out.println("Quantidade de Eigenfaces utilizadas:  " + vectorNumbers);

        for (int img = 0; img < linearizedImages[0].length; img++) {
            double[] columnWeights = getColumn(weightsMatrix, img);
            double[] reconstructedImage = reconstructImage(averageVectors, eigenfaces, columnWeights, vectorNumbers);
            double[][] reconstructedImageMatrix = array1DToMatrix(reconstructedImage, allMatricesCsv[img]);
            double maximumAbsolutError = calculateMAE(allMatricesCsv[img], reconstructedImageMatrix);
            System.out.println("Para a imagem: " + csvFiles[img] + ", foi utilizado este vetor peso : " + Arrays.toString(columnWeights));
            System.out.println("O erro absoluto médio dessa imagem com sua original foi: %.3f" + maximumAbsolutError);
            saveImage(reconstructedImageMatrix, csvFiles[img], "Output/Func2/ImagensReconstruidas", 0);
            saveMatrixToFile(reconstructedImageMatrix, csvFiles[img], "Output/Func2/Eigenfaces", 0);
        }
    }

    public static void calculateFunction3(int vectorNumbers, String csvLocation, String dataBase) {
        String[] csvFiles = getCSVFileNames(dataBase);
        double[][][] allMatricesCsv = getMatricesFromCsvFolder(dataBase);
        double[][] oneMatrixCsv = readCSVToMatrix(csvLocation);


        double[][] linearizedImages = new double[allMatricesCsv[0].length * allMatricesCsv[0].length][allMatricesCsv.length];
        double[][] weightsMatrix = new double[allMatricesCsv[0].length * allMatricesCsv[0].length][allMatricesCsv.length];
        populateLinearizedImages(linearizedImages, allMatricesCsv);
        double[] averageVectors = calculateMeanVector(linearizedImages);
        double[][] phi = centralizeImages(linearizedImages, averageVectors);
        int vectorK = validateEigenVectors(linearizedImages, vectorNumbers);

        double[][] phiT = transposeMatrix(phi);
        double[][] phiTxPhi = multiplyMatrices(phiT, phi);
        double[][] eigenVectors = getEigenVectors(phiTxPhi);
        double[][] selectedColumnsK = getValuesAndIndexArray(eigenVectors, vectorK);
        double[][] newEigenVectorsK = createSubmatrix(eigenVectors, selectedColumnsK);
        double[][] expandedVectorsK = multiplyMatrices(phi, newEigenVectorsK);
        double[][] eigenfaces = normalize(expandedVectorsK);

        populateWeightsMatrix(weightsMatrix, phi, eigenfaces);

        identifyClosestImage(vectorK, csvFiles, averageVectors, eigenfaces, oneMatrixCsv, weightsMatrix, allMatricesCsv);
    }

    public static void identifyClosestImage(int vectorNumbers, String[] csvFiles, double[] averageVectors, double[][] eigenfaces, double[][] oneMatrixCsv, double[][] weightsMatrix, double[][][] allMatricesCsv) {
        double[] linearizedPrincipalImage = matrixToArray1D(oneMatrixCsv);
        double[] phiPrincipalImage = subtractionColumns(linearizedPrincipalImage, averageVectors);

        double[] principalWeightsVector = calculateWeights(phiPrincipalImage, eigenfaces);

        double[] distances = calculateEuclidianDistance(principalWeightsVector, weightsMatrix);
        int closestImageIndex = checkCloserVetor(distances);

        double[] closestImageWeights = getColumn(weightsMatrix, closestImageIndex);
        double[] reconstructedImage = reconstructImage(averageVectors, eigenfaces, closestImageWeights, vectorNumbers);

        double[][] reconstructedImageMatrix = array1DToMatrix(reconstructedImage, allMatricesCsv[0]);

        printFunction3(vectorNumbers, csvFiles, closestImageIndex, distances);
        saveImage(reconstructedImageMatrix, csvFiles[closestImageIndex], "Output/Func3/Identificacao", 1);
    }

    public static void generateNewImage(int vectorNumbers, String dataBase) {
        double[][][] allMatricesCsv = getMatricesFromCsvFolder(dataBase);
        double[][] linearizedImages = new double[allMatricesCsv[0].length * allMatricesCsv[0].length][allMatricesCsv.length];
        double[][] weightsMatrix = new double[allMatricesCsv[0].length * allMatricesCsv[0].length][allMatricesCsv.length];
        populateLinearizedImages(linearizedImages, allMatricesCsv);
        double[] meanVector = calculateMeanVector(linearizedImages);
        double[][] phi = centralizeImages(linearizedImages, meanVector);
        int vectorK = validateEigenVectors(linearizedImages, vectorNumbers);

        double[][] phiT = transposeMatrix(phi);
        double[][] phiTxPhi = multiplyMatrices(phiT, phi);
        double[][] eigenVectors = getEigenVectors(phiTxPhi);
        double[][] eigenValues = getEigenValues(phiTxPhi);
        double[][] selectedColumnsK = getValuesAndIndexArray(eigenValues, vectorK);
        double[][] newEigenVectorsK = createSubmatrix(eigenVectors, selectedColumnsK);
        double[][] newEigenValuesK = constructDiagonalMatrix(selectedColumnsK);
        double[][] expandedVectorsK = multiplyMatrices(phi, newEigenVectorsK);
        double[][] eigenfaces = normalize(expandedVectorsK);

        populateWeightsMatrix(weightsMatrix, phi, eigenfaces);
        int dimension = meanVector.length;
        double[] newImage = creationImage(dimension, meanVector, vectorK, newEigenValuesK, eigenfaces);
        double[][] newImageMatrix = array1DToMatrix(newImage, allMatricesCsv[0]);
        saveImage(newImageMatrix, "Input/Funcao2-3/csv", "Output/Func4", 0);
    }
    //* ------------------ Fim dos métodos de distribuição de tarefas ------------------


    //* ------------------ Métodos de cálculos ------------------
    public static EigenDecomposition decomposeMatrix(double[][] matrixToDecompose) {
        Array2DRowRealMatrix decomposedMatrix = new Array2DRowRealMatrix(matrixToDecompose);
        return new EigenDecomposition(decomposedMatrix);
    }

    public static void quitApplication() {
        uiQuitParameterMenu();
        receiveExitConfirmation(null);
    }

    public static void populateWeightsMatrix(double[][] weightsMatrix, double[][] phi, double[][] eigenfaces) {
        for (int img = 0; img < phi[0].length; img++) {
            double[] actualPhiColumn = getColumn(phi, img);
            double[] weights = calculateWeights(actualPhiColumn, eigenfaces);

            for (int i = 0; i < weights.length; i++) {
                weightsMatrix[i][img] = weights[i];
            }
        }
    }

    public static void populateLinearizedImages(double[][] linearizedImages, double[][][] imageMatrices) {
        for (int img = 0; img < imageMatrices.length; img++) {
            double[] linearizedMatrix = matrixToArray1D(imageMatrices[img]);
            for (int i = 0; i < linearizedMatrix.length; i++) {
                linearizedImages[i][img] = linearizedMatrix[i];
            }
        }
    }

    public static int validateEigenVectors(double[][] matrix, int vectorNumbers) {
        if (vectorNumbers <= MIN_QUANTITY_VECTORS || vectorNumbers > matrix[0].length) {
            vectorNumbers = matrix[0].length;
        }
        return vectorNumbers;
    }

    public static double calculateMAE(double[][] originalMatrix, double[][] matrixEigenFaces) {
        int rows = originalMatrix.length;
        int columns = originalMatrix[0].length;
        double errorAbsMed = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                errorAbsMed += Math.abs(originalMatrix[i][j] - matrixEigenFaces[i][j]);
            }
        }
        return errorAbsMed / (rows * columns);
    }

    public static double[] getColumn(double[][] matrix, int column) {
        double[] columnData = new double[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            columnData[i] = matrix[i][column];
        }
        return columnData;
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

    public static double[] calculateWeights(double[] phi, double[][] eigenfaces) {
        if (phi.length != eigenfaces.length) {
            errorGeneral("Para calcular os pesos o comprimento de 'phi' deve ser igual a quantidade de linhas da matriz 'eigenfaces'.");
        }

        double[] weights = new double[phi.length];

        for (int j = 0; j < eigenfaces[0].length; j++) {
            weights[j] = 0;
            for (int i = 0; i < eigenfaces.length; i++) {
                weights[j] += phi[i] * eigenfaces[i][j];
            }
        }
        return weights;
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

    public static double[][] centralizeImages(double[][] images, double[] meanVector) {
        int numPixels = meanVector.length;
        int numImages = images[0].length;
        if (images.length != numPixels) {
            errorGeneral("Para centralizar a imagem o número de pixels na matriz de imagens deve ser igual ao tamanho do vetor médio.");
        }

        double[][] phi = new double[numPixels][numImages];

        for (int i = 0; i < numPixels; i++) {
            for (int j = 0; j < numImages; j++) {
                phi[i][j] = images[i][j] - meanVector[i];
            }
        }

        return phi;
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

    public static double[][] constructDiagonalMatrix(double[][] matrixvaluesK) {
        double[][] matrixvaluesKPrint = new double[matrixvaluesK.length][matrixvaluesK.length];
        for (int i = 0; i < matrixvaluesK.length; i++) {
            matrixvaluesKPrint[i][i] = matrixvaluesK[i][0];
        }
        return matrixvaluesKPrint;
    }

    public static double[][] getValuesAndIndexArray(double[][] eigenValuesArray, int eigenfaces) {
        double[][] valuesAndIndexArray = new double[eigenfaces][2];

        for (int i = 0; i < valuesAndIndexArray.length; i++) {
            valuesAndIndexArray[i][0] = Double.MIN_VALUE;
        }

        for (int i = 0; i < eigenValuesArray.length; i++) {
            double absValue = Math.abs(eigenValuesArray[i][i]);
            for (int j = 0; j < valuesAndIndexArray.length; j++) {
                if (absValue > Math.abs(valuesAndIndexArray[j][0])) {
                    // Joga os valores para a direita se encontrar um valor maior mais ao fim da matrix
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

    public static double[][] array1DToMatrix(double[] reconstructedImage, double[][] matrixSizeExample) {
        int rows = matrixSizeExample.length;
        int columns = matrixSizeExample[0].length;
        if (reconstructedImage.length != rows * columns) {
            errorGeneral("O tamanho do vetor não corresponde às dimensões da matriz.");
        }

        double[][] reconstructedMatrix = new double[rows][columns];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                reconstructedMatrix[i][j] = reconstructedImage[i * columns + j];
            }
        }
        return reconstructedMatrix;
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

    public static double[] reconstructImage(double[] averageVector, double[][] eigenfaces, double[] columnWeights, int quantityEigenfaces) {
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

    public static double[] calculateEuclidianDistance(double[] principalVector, double[][] weightsMatrix) {
        if (principalVector.length != weightsMatrix.length) {
            throw new IllegalArgumentException("O comprimento do vetor principal não corresponde ao número de linhas da matriz de pesos.");
        }

        double[] result = new double[weightsMatrix[0].length];
        for (int i = 0; i < weightsMatrix[0].length; i++) {
            double sum = 0;
            for (int j = 0; j < weightsMatrix.length; j++) {
                sum += Math.pow(principalVector[j] - weightsMatrix[j][i], 2);
            }
            result[i] = Math.sqrt(sum);
        }
        return result;
    }

    public static int checkCloserVetor(double[] distances) {
        double minDistance = Double.MAX_VALUE;
        int closestImageIndex = -1;

        for (int i = 0; i < distances.length; i++) {
            if (distances[i] < minDistance) {
                minDistance = distances[i];
                closestImageIndex = i;
            }
        }
        return closestImageIndex;
    }

    public static double[] creationImage(int dimension, double[] meanVector,int k, double[][] lambdas, double[][] eigenfaces) {
        double[] newImage = new double[dimension];
        for (int i = 0; i < dimension; i++) {
            newImage[i] = meanVector[i];
        }
        for (int i = 0; i < k; i++) {
            double weightImage = Math.random() * (2 * Math.sqrt(lambdas[i][i])) - Math.sqrt(lambdas[i][i]);
            for (int j = 0; j < dimension; j++) {
                newImage[j] += weightImage * eigenfaces[i][j];
            }
        }

        return newImage;
    }
    //* ------------------ Fim dos métodos de cálculos ------------------


    //* ------------------ Métodos de entrada e saída ------------------
    public static void writeArrayAsImage(int[][] array, String outputFilePath) throws IOException {
        int height = array.length;
        int width = array[0].length;

        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int intensity = array[y][x];
                if (intensity < MIN_BIT_VALUE || intensity > MAX_BIT_VALUE) {
                    errorGeneral("Erro: Na normalização dos pixels, a intensidade do pixel deve estar entre 0 e 255.");
                }
                int rgb = (intensity << 16) | (intensity << 8) | intensity;
                image.setRGB(x, y, rgb);
            }
        }

        File outputFile = new File(outputFilePath);
        ImageIO.write(image, "jpg", outputFile);
    }

    public static void saveImage(double[][] imageArray, String inputCsvPath, String outputFolderPath, int printOrNot) {
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

        String jpgFileName = new File(inputCsvPath).getName();
        if (inputCsvPath.endsWith(".csv")) {
            // Substituir extensão .csv por .jpg
            jpgFileName = jpgFileName.replace(".csv", ".jpg");
        } else {
            // Se for uma pasta, usar o nome da pasta + "_output.jpg"
            String folderName = inputCsvPath.substring(inputCsvPath.lastIndexOf("/") + 1);
            jpgFileName = folderName + "_output.jpg";
        }
        String outputPath = outputFolderPath + "/" + jpgFileName;


        File outputFolder = new File(outputFolderPath);
        if (!outputFolder.exists()) {
            if (!outputFolder.mkdirs()) {
                System.err.println("Falha ao criar o diretório: " + outputFolderPath);
                return;
            }
        }

        try {
            writeArrayAsImage(normalizedImage, outputPath);
            if (printOrNot == 1) {
                System.out.println("A imagem foi gerada com sucesso: " + outputPath);
            }
        } catch (IOException e) {
            System.err.println("Erro ao salvar a imagem: " + e.getMessage());
        }
    }

    public static void saveMatrixToFile(double[][] matrix, String inputCsvPath, String outputFolderPath, int printOrNot) {
        String csvFileName = new File(inputCsvPath).getName();
        String newFileName = "Reconstruct-" + csvFileName;

        File file = new File(outputFolderPath + "/" + newFileName);


        try (PrintWriter writer = new PrintWriter(file)) {
            for (double[] row : matrix) {
                String rowString = String.join(" , ", Arrays.stream(row)
                        .mapToObj(val -> String.format("%.0f", val))
                        .toArray(String[]::new));
                writer.println(rowString);
            }
            if (printOrNot == 1) {
                System.out.println("Arquivo CSV criado com sucesso: " + file.getName());
            }
        } catch (IOException e) {
            System.err.println("Erro ao salvar a matriz no arquivo: " + e.getMessage());
            e.printStackTrace();
        }
    }

    public static String receiveDataBaseLocation(String[] args) {
        String dataBaseLocationArgs;
        if (args == null) {
            String dataBaseLocation = scanner.next();
            if (!checkDataBaseLocation(dataBaseLocation)) {
                System.out.println("Erro: Localização da base de imagens inválida.");
                System.out.println("Tentar novamente ? (S/N)");
                String answer;
                do {
                    answer = scanner.next().toUpperCase();
                    if (answer.equals("S")) {
                        dataBaseLocation = verifyDataBaseLocation();
                    } else if (answer.equals("N")) {
                        System.out.println("Saindo da aplicação, ainda poderá retornar ao menu inicial.");
                        quitApplication();
                    } else {
                        System.out.println("Opção inválida, responda com S/N.");
                    }
                } while (!answer.equals("S") && !answer.equals("N"));
            }
            return dataBaseLocation;
        } else {
            dataBaseLocationArgs = args[5];
            if (!checkDataBaseLocation(dataBaseLocationArgs)) {
                errorGeneral("Erro: Localização da base de dados inválida.");
            }
            return dataBaseLocationArgs;
        }
    }

    public static String receiveCsvLocation(String[] args) {
        String csvLocationArgs;
        if (args == null) {
            String csvLocation = scanner.next();
            if (!checkCsvLocation(csvLocation)) {
                System.out.println("Erro: Localização do csv inválida");
                System.out.println("Tentar novamente ? (S/N)");
                String answer;
                do {
                    answer = scanner.next().toUpperCase();
                    if (answer.equals("S")) {
                        csvLocation = verifyCsvLocation();
                    } else if (answer.equals("N")) {
                        System.out.println("Saindo da aplicação, ainda poderá retornar ao menu inicial.");
                        quitApplication();
                    } else {
                        System.out.println("Opção inválida, responda com S/N.");
                    }
                } while (!answer.equals("S") && !answer.equals("N"));
            }
            return csvLocation;
        } else {
            csvLocationArgs = args[5];
            if (!checkCsvLocation(csvLocationArgs)) {
                errorGeneral("Erro: Localização do csv inválida");
            }
            return csvLocationArgs;
        }
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

    public static double[][] readCSVToMatrix(String path) {
        try {
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

            double[][] matrix = new double[rowCount][columnCount];

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
    //* ------------------ Fim dos Métodos de entrada e saída ------------------


    //* ------------------ Verificações ------------------
    public static boolean checkCorrectParametersStructure(String[] parameters) {
        if (parameters.length == 8) {
            return parameters[0].equals("-f") && parameters[2].equals("-k") && parameters[4].equals("-i") && parameters[6].equals("-j");
        }
        return false;
    }

    public static boolean checkFunctionOptions(int function) {
        return function >= 1 && function <= 5 || function == 0;
    }

    public static boolean checkSizeBoundaries(int rows, int cols) {
        return rows > MAX_SIZE_ROWS_AND_COLS || cols > MAX_SIZE_ROWS_AND_COLS || rows < MIN_SIZE_ROWS_AND_COLS || cols < MIN_SIZE_ROWS_AND_COLS;
    }

    public static boolean checkCsvLocation(String csvLocation) {
        File csv = new File(csvLocation);
        if (csvLocation.isEmpty()) {
            return false;
        } else if (!csvLocation.contains(".csv")) {
            return false;
        } else return csv.exists();
    }

    public static boolean checkDataBaseLocation(String dataBaseLocation) {
        File imageDirectory = new File(dataBaseLocation);
        if (dataBaseLocation.isEmpty()) {
            return false;
        }
        return imageDirectory.exists();
    }

    public static String verifyCsvLocation() {
        String csvLocation;
        uiCsvLocation();
        csvLocation = receiveCsvLocation(null);
        return csvLocation;
    }

    public static String verifyDataBaseLocation() {
        String dataBaseLocation;
        do {
            uiDataBase();
            dataBaseLocation = receiveDataBaseLocation(null);
        } while (!checkDataBaseLocation(dataBaseLocation));
        return dataBaseLocation;
    }

    public static int verifyFunction() {
        int function;
        do {
            uiInitialMenu();
            function = receiveFunction(null);
        } while (!checkFunctionOptions(function));
        return function;
    }

    public static int verifyVectorNumbers() {
        int vectorNumbers;
        do {
            uiVectorNumbers();
            vectorNumbers = receiveNumberVectors(null);
        } while (vectorNumbers < MIN_QUANTITY_VECTORS);
        return vectorNumbers;
    }

    public static void checkExistanceFileDirectory(String csvLocation) {
        try {
            scanner = new Scanner(new File(csvLocation));
        } catch (FileNotFoundException e) {
            errorGeneral("Erro ao abrir os arquivos: " + e.getMessage());
        }
    }
    //* ------------------ Fim verificações ------------------


    //* ------------------ Métodos de Interfaces ------------------
    public static void uiInitialMenu() {
        System.out.print("\n+----------------------------------------------------+\n");
        System.out.println("|            Qual função deseja realizar?            |");
        System.out.println("+----------------------------------------------------+");
        System.out.println("|                                                    |");
        System.out.println("| 1 - Decomposição própria de uma matriz simétrica.  |");
        System.out.println("| 2 - Reconstrução de imagens usando Eigenfaces.     |");
        System.out.println("| 3 - Identificação de imagem mais próximas.         |");
        System.out.println("| 4 - Gerar uma imagem aleatória com Eigenfaces.     |");
        System.out.println("| 5 - Conheça a equipa de desenvolvimento!           |");
        System.out.println("| 0 - Deseja sair da aplicação ?                     |");
        System.out.println("+----------------------------------------------------+");
        System.out.print("Opção: ");
    }

    public static void uiDevTeam() {
        System.out.println("+----------------------------------------------------+");
        System.out.println("|             Equipa de desenvolvimento:             |");
        System.out.println("+----------------------------------------------------+");
        System.out.println("|Alexandre Pereira Henrique                          |");
        System.out.println("|Luiz Gabriel de Souza Sargaço Teixeira              |");
        System.out.println("|Rafael Pinto Vieira                                 |");
        System.out.println("|Rita Mafalda Martins de Oliveira                    |");
        System.out.println("|Younés André Marques de Almeida Bouayad             |");
        System.out.println("+----------------------------------------------------+");
    }

    public static void uiVectorNumbers() {
        System.out.println("----- Quantos vetores próprios deseja utilizar? -----");
        System.out.print("Quantidade: ");
    }

    public static void uiCsvLocation() {
        System.out.println("-- Qual a localização do csv que deseja utilizar? --");
        System.out.print("Localização: ");
    }

    public static void uiDataBase() {
        System.out.println("------ Qual a localização da base de imagens? ------");
        System.out.print("Localização: ");
    }

    public static void uiQuitParameterMenu() {
        System.out.println("-- Tem certeza que deseja sair da aplicação? (S/N) --");
        System.out.print("Opção: ");
    }
    //* --------------------- Fim menus de opções ------------------


    //* ------------------ Receber parâmetros ------------------
    public static int receiveFunction(String[] args) {
        int functionArgs;
        if (args == null) {
            int function = scanner.nextInt();
            if (!checkFunctionOptions(function)) {
                System.out.println("Erro: Opção inválida.");
                System.out.println("Tente novamente.");
            }
            return function;
        } else {
            functionArgs = Integer.parseInt(args[1]);
            if (!checkFunctionOptions(functionArgs)) {
                errorGeneral("Erro: Opção inválida, as opções são: 1 a 4, opções como 5 e 0 estão disponíveis apenas para o modo interativo!.");
            }
            return functionArgs;
        }
    }

    public static int receiveNumberVectors(String[] args) {
        int vectorNumbersArgs;
        if (args == null) {
            vectorNumbersArgs = scanner.nextInt();
            if (vectorNumbersArgs < MIN_QUANTITY_VECTORS) {
                System.out.println("Erro: O número de vetores deve ser maior que 0.");
                System.out.println("Tentar novamente? (S/N)");
                String answer = scanner.next().toUpperCase();
                if (answer.equals("N")) {
                    System.out.println("A sair da aplicação, ainda poderá desistir mas\nretornará ao menu inicial.");
                    quitApplication();
                }
            }
        } else {
            vectorNumbersArgs = Integer.parseInt(args[3]);
            if (vectorNumbersArgs < MIN_QUANTITY_VECTORS) {
                errorGeneral("Erro: O número de vetores deve ser maior que 0.");
            }
            return vectorNumbersArgs;
        }
        return vectorNumbersArgs;
    }

    public static void receiveExitConfirmation(String[] args) {
        if (args == null) {
            String confirmeExit;
            do {
                confirmeExit = scanner.next().toUpperCase();
                if (confirmeExit.equals("S")) {
                    System.exit(0);
                } else if (confirmeExit.equals("N")) {
                    System.out.println("A retornar para o menu inicial.");
                    System.out.println();
                    runInterative();
                } else {
                    System.out.println("Erro: Responda com S/N");
                }
            } while (!confirmeExit.equals("N"));
        }
    }
    //* -------------------- Fim receber parâmetros ------------------


    //* ------------------ Operações básicas com Matrizes ------------------
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

    public static double[][] transposeMatrix(double[][] matrix) {
        double[][] transposedMatrix = new double[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                transposedMatrix[j][i] = matrix[i][j];
            }
        }
        return transposedMatrix;
    }

    public static double[] subtractionColumns(double[] columnLeft, double[] columnRight) {
        double[] matrixResult = new double[columnLeft.length];
        for (int i = 0; i < columnLeft.length; i++) {
            matrixResult[i] = columnLeft[i] - columnRight[i];
        }

        return matrixResult;
    }

    public static double[][] createSubmatrix(double[][] eigenVectors, double[][] valuesAndIndexArray) {
        boolean[] keepColumnsBoolean = new boolean[eigenVectors[0].length];

        for (double[] columns : valuesAndIndexArray) {
            keepColumnsBoolean[(int) columns[1]] = true;
        }

        double[][] submatrix = new double[eigenVectors.length][valuesAndIndexArray.length];

        int subMatrixRows = 0;
        for (double[] doubles : eigenVectors) {
            int subMatrixColumns = 0;
            for (int j = 0; j < doubles.length; j++) {
                if (keepColumnsBoolean[j]) {
                    submatrix[subMatrixRows][subMatrixColumns] = doubles[j];
                    subMatrixColumns++;
                }
            }
            subMatrixRows++;
        }
        return submatrix;
    }
    //* ----------------- Fim operações básicas com matrizes ------------------


    //* -------------------- Printar Matrizes -----------------------
    public static void printMatrix(double[][] matrixToPrint, String matrixName) {
        System.out.println("\nMatriz: " + matrixName + " ↓");
        printLine(matrixToPrint[0].length, "____________");
        System.out.println();

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

    public static void printLine(int length, String pattern) {
        for (int i = 0; i < length; i++) {
            System.out.print(pattern);
        }
    }
    //* ----------------- Fim Printar matrizes -----------------------


    //* ----------------- Printar Funcionalidades -----------------
    public static void printFunction1(int numbersEigenfaces, double[][] newEigenValuesK, double[][] newEigenVectorsK, double maximumAbsolutError) {
        System.out.println("A quantidade de Eigenfaces selecionadas para a variável K foi: " + numbersEigenfaces);
        printMatrix(newEigenValuesK, "Valores Próprios da matriz K");
        printMatrix(newEigenVectorsK, "Vetores Próprios matriz K:");
        System.out.printf("Erro Absoluto Médio: %.3f\n", maximumAbsolutError);
    }

    public static void printFunction3(int vectorNumbers, String[] csvFiles, int closestImageIndex, double[] distances) {
        System.out.println("O número de vetores próprios utilizados: " + vectorNumbers);
        System.out.printf("A imagem mais próxima foi: %s e foi salva em Identificação!\n\n", csvFiles[closestImageIndex]);

        for (int i = 0; i < csvFiles.length; i++) {
            if (i == closestImageIndex) {
                System.out.printf("Essa foi a imagem mais próxima da solicitada! %s e sua distância foi: %.1f\n", csvFiles[i], distances[i]);
            } else {
                System.out.printf("Distância euclidiana para a imagem %s: %.1f\n", csvFiles[i], distances[i]);
            }
        }
    }

    public static void printHeaderFunction(String functionName) {
        int length = functionName.length();
        int totalWidth = length + 4; // 2 espaços de cada lado do texto
        int padding = (totalWidth - length) / 2;

        System.out.print("\n+");
        printLine(totalWidth, "-");
        System.out.print("+\n");
        System.out.printf("|%" + padding + "s%s%" + padding + "s|\n+", "", functionName, "");
        printLine(totalWidth, "-");
        System.out.print("+\n");
        System.out.println();
    }
    //* ----------------- Fim printar funcionalidades -----------------------


    //! ------------------ Error Messages para não interativo ------------------
    public static void errorGeneral(String error) {
        //! Esse tipo de mensagem de erro deve ser usado apenas para o modo não interativo!
        System.out.println(error);
        System.exit(1);
    }
    //! ------------------ Fim error messages --------------
}