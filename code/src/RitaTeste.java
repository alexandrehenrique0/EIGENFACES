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
    // constantes para limites de tamanho.
    public static final int MAX_SIZE_ROWS = 256;
    public static final int MAX_SIZE_COLS = 256;
    public static final int MIN_SIZE_ROWS = 1;
    public static final int MIN_SIZE_COLS = 1;
    public static final int MIN_QUANTITY_VECTORS = 1;

    // scanners para inputs.
    public static Scanner scanner = new Scanner(System.in);
    public static Scanner scannerCsv = new Scanner(System.in);

    public static void main(String[] args) {
        // Verificação de parâmetros para decidir entre interativo e não interativo.
        if (checkCorrectParametersStructure(args)) {
            runNonInteractive(args);
        } else {
            runInterative();
        }
        scanner.close();
        scannerCsv.close();
    }

    //* ------------------ Modos de execução ------------------
    public static void runInterative() {
        // Parametros de entrada
        int function = 0;
        int vectorNumbers = 0;
        String csvLocation;
        String imageFolderLocation;

        // Roda enquanto a função for inválida
        function = verifyFunction();

        // Para sair da aplicação
        if (function == 4) {
            quitApplication();
        }

        vectorNumbers = verifyVectorNumbers();
        csvLocation = verifyCsvLocation();
        imageFolderLocation = verifyImageFolderLocation(function);

        //? esse metodo checkExistanceFileDirectory(csvLocation, imageFolderLocation) já
        //? não é testado em verifyCsvLocation() e verifyImageFolderLocation() ?

        //Obtém os nomes dos arquivos CSV em um diretório de imagens e retorna uma lista de aqrquivos CSV localizados no diretório especifico
        String[] csvFiles = getCSVFileNames(imageFolderLocation);

        //Lê um arquivo CSV no caminho fornecido e converte seu conteúdo em uma matriz 2D (double[][])
        double[][] oneMatrixCsv = readCSVToMatrix(csvLocation);

        //retorna várias matrizes 2D de diferentes arquivos CSV dentro do diretório.
        double[][][] allMatricesCsv = getMatricesFromCsvFolder(imageFolderLocation);

        switchPrimaryFunctions(function, vectorNumbers, csvLocation, csvFiles, oneMatrixCsv, allMatricesCsv);
    }

    public static void runNonInteractive(String[] args) {
        // Parametros de entrada
        int function;
        int vectorNumbers;
        String csvLocation;
        String imageFolderLocation;

        // Receber os parâmetros
        function = receiveFunction(args);
        vectorNumbers = receiveNumberVectors(args);
        csvLocation = receiveCsvLocation(args);
        imageFolderLocation = receiveImageLocation(args);

        // Verificar se os arquivo e diretório existem
        checkExistanceFileDirectory(csvLocation);

        String[] csvFiles = getCSVFileNames(imageFolderLocation);

        // Obter a matriz do CSV
        double[][] oneMatrixCsv = readCSVToMatrix(csvLocation);

        // Obter a matriz do CSV para a função 2
        double[][][] allMatricesCsv = getMatricesFromCsvFolder(imageFolderLocation);

        // Função que contém as funções principais
        switchPrimaryFunctions(function, vectorNumbers, csvLocation, csvFiles, oneMatrixCsv, allMatricesCsv);
    }
    //* ------------------ Fim modos de execução ------------------


    //* ------------------ Métodos principais ------------------
    public static void switchPrimaryFunctions(int function, int vectorNumbers, String csvLocation, String[] csvFiles, double[][] oneMatrixCsv, double[][][] allMatricesCsv) {
        // Common variables for all 3 functions
        double[][] linearizedImages = new double[allMatricesCsv[0].length * allMatricesCsv[0].length][allMatricesCsv.length];
        double[][] weightsMatrix = new double[allMatricesCsv[0].length * allMatricesCsv[0].length][allMatricesCsv.length];
        populateLinearizedImages(linearizedImages, allMatricesCsv);
        double[] averageVectors = calculateAverageVector(linearizedImages);
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
                decomposeSymmetricMatrix(oneMatrixCsv, vectorNumbers, csvLocation);
                System.out.println();
                System.out.println("Funcionalidade 1 finalizada.");
                runInterative();
                break;
            case 2:
                printHeaderFunction("Reconstrução de Imagens usando Eigenfaces");
                reconstructImagesWithEigenfaces(vectorK, csvFiles, averageVectors, eigenfaces, linearizedImages, weightsMatrix, allMatricesCsv);
                System.out.println();
                System.out.println("Funcionalidade 2 finalizada.");
                runInterative();
                break;
            case 3:
                printHeaderFunction("Identificação de imagem mais próxima");
                identifyClosestImage(vectorK, csvFiles, averageVectors, oneMatrixCsv, weightsMatrix, allMatricesCsv, eigenfaces);
                System.out.println();
                System.out.println("Funcionalidade 3 finalizada.");
                runInterative();
                break;
        }
    }

    public static void decomposeSymmetricMatrix(double[][] oneMatrixCsv, int vectorNumbers, String csvLocation) {

        double[][] eigenVectors = getEigenVectors(oneMatrixCsv);
        double[][] eigenValues = getEigenValues(oneMatrixCsv);

        if (oneMatrixCsv[0].length < vectorNumbers) {
            vectorNumbers = oneMatrixCsv[0].length;
        }

        double[][] valuesAndIndexArray = getValuesAndIndexArray(eigenValues, vectorNumbers);
        double[][] newEigenVectorsK = createSubmatrix(eigenVectors, valuesAndIndexArray);
        double[][] newEigenValuesK = constructDiagonalMatrix(valuesAndIndexArray);
        double[][] newEigenVectorsTransposeK = transposeMatrix(newEigenVectorsK);
        double[][] matrixEigenFaces = multiplyMatrices(multiplyMatrices(newEigenVectorsK, newEigenValuesK), newEigenVectorsTransposeK);

        double maximumAbsolutError = calculateMAE(oneMatrixCsv, matrixEigenFaces);

        saveMatrixToFile(matrixEigenFaces, csvLocation, "Output/Func1", 1);
        printFunction1(vectorNumbers, newEigenValuesK, newEigenVectorsK, maximumAbsolutError);
    }

    public static void reconstructImagesWithEigenfaces(int vectorNumbers, String[] csvFiles, double[] averageVectors, double[][] eigenfaces, double[][] linearizedImages, double[][] weightsMatrix, double[][][] allMatricesCsv) {

        System.out.println("Valores do vetor médio: " + Arrays.toString(averageVectors));
        System.out.println("Quantidade de Eigenfaces utilizadas:  " + vectorNumbers);

        for (int img = 0; img < linearizedImages[0].length; img++) {
            double[] columnWeights = getColumn(weightsMatrix, img);
            double[] reconstructedImage = reconstructImage(averageVectors, eigenfaces, columnWeights, vectorNumbers);
            double[][] reconstructedImageMatrix = array1DToMatrix(reconstructedImage, allMatricesCsv[img]);
            System.out.println("Para a imagem: " + csvFiles[img] + ", foi utilizado este vetor peso : " + Arrays.toString(columnWeights));
            saveImage(reconstructedImageMatrix, csvFiles[img], "Output/Func2/ImagensReconstruidas", 0);
            saveMatrixToFile(reconstructedImageMatrix, csvFiles[img], "Output/Func2/Eigenfaces", 0);
        }
    }

    public static void identifyClosestImage(int vectorNumbers, String[] csvFiles, double[] averageVectors, double[][] oneMatrixCsv, double[][] weightsMatrix, double[][][] allMatricesCsv, double[][] eigenfaces) {
        double[] linearizedPrincipalImage = matrixToArray1D(oneMatrixCsv);
        double[] phiPrincipalImage = subtractionColumns(linearizedPrincipalImage, averageVectors);

        double[] principalWeightsVector = calculateWeights(phiPrincipalImage, eigenfaces);

        double[] distances = calculateEuclidianDistance(principalWeightsVector, weightsMatrix);
        int closestImageIndex = checkCloserVetor(distances);

        double[] closestImageWeights = getColumn(weightsMatrix,closestImageIndex);
        double[] reconstructedImage = reconstructImage(averageVectors, eigenfaces, closestImageWeights, vectorNumbers);

        System.out.println("O número de vetores próprios utilizados: " + vectorNumbers);
        System.out.printf("A imagem mais próxima foi: %s e foi salva em Identificação!\n\n", csvFiles[closestImageIndex]);

        double[][] reconstructedImageMatrix = array1DToMatrix(reconstructedImage, allMatricesCsv[0]);

        for (int i = 0; i < csvFiles.length; i++) {
            if (i == closestImageIndex) {
                System.out.printf("Essa foi a imagem mais próxima da solicitada! %s e sua distância foi: %.1f\n", csvFiles[i], distances[i]);
            } else {
                System.out.printf("Distância euclidiana para a imagem %s: %.1f\n", csvFiles[i], distances[i]);
            }
        }

        saveImage(reconstructedImageMatrix, csvFiles[closestImageIndex], "Output/Func3/Identificacao", 1);
    }
    //* ------------------ fim dos metodos principais ------------------


    //* ------------------ Funcionalidades comuns ------------------
    public static EigenDecomposition decomposeMatrix(double[][] matrixToDecompose) {
        Array2DRowRealMatrix decomposedMatrix = new Array2DRowRealMatrix(matrixToDecompose);
        return new EigenDecomposition(decomposedMatrix);
    }

    public static void quitApplication() {
        uiQuitParameterMenu();
        receiveExitConfirmation(null);
    }

    public static void populateRow(double[][] matrix, int row, String line) {
        String[] values = line.split(",");
        for (int col = 0; col < values.length; col++) {
            try {
                matrix[row][col] = Double.parseDouble(values[col].trim());
            } catch (NumberFormatException e) {
                matrix[row][col] = 0; // or any default value you prefer
            }
        }
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

    public static int[] getDimensions() {
        int rows = 0;
        int cols = 0;
        while (scannerCsv.hasNextLine()) {
            String line = scannerCsv.nextLine().trim();
            if (!line.isEmpty()) {
                if (rows == 0) {
                    cols = line.split(",").length;
                }
                rows++;
            }
        }
        if (checkSizeBoundaries(rows, cols)) {
            errorGeneral("Erro: Dimensões da matriz fora dos limites: " + rows + "x" + cols);
        }

        return new int[]{rows, cols};
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
            errorGeneral("O comprimento de 'phi' deve ser igual ao número de linhas em 'eigenfaces'.");
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

    public static double[][] centralizeImages(double[][] images, double[] meanVector) {
        int numPixels = meanVector.length;
        int numImages = images[0].length;
        if (images.length != numPixels) {
            errorGeneral("O número de pixels na matriz de imagens deve ser igual ao tamanho do vetor médio.");
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

    public static double[][] covariances(double[][] matrixA) {
        int quantityOfImages = matrixA[0].length;
        double[][] matrixATransposed = transposeMatrix(matrixA);
        double[][] matrixATmultiplyByA = multiplyMatrices(matrixATransposed, matrixA);
        return multiplyMatrixEscalar(matrixATmultiplyByA, 1.0 / quantityOfImages);
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
    //* ------------------ Fim funcionalidades comuns ------------------


    //* -------------------- Exclusivo funcionalidade 1 -----------------------
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
    //* ----------------- Fim funcionalidade 1 ------------------


    //* ----------------- Exclusivo funcionalidade 2 ------------------
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

    public static void populateMatrix(double[][] matrix, String csvLocation) {
        try {
            scannerCsv = new Scanner(new File(csvLocation));
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Erro ao reabrir o arquivo CSV: " + e.getMessage());
        }
        int row = 0;
        while (scannerCsv.hasNextLine()) {
            String line = scannerCsv.nextLine().trim();
            if (!line.isEmpty()) {
                populateRow(matrix, row, line);
                row++;
            }
        }
        scannerCsv.close();
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
    //* ----------------- Fim das funcionalidades 2 ------------------


    //* ----------------- Exclusivo funcionalidade 3 ------------------
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
    //* ----------------- Fim funcionalidade 3 ------------------


    //* ------------------ Metodos de entrada e saída ------------------
    public static void writeArrayAsImage(int[][] array, String outputFilePath) throws IOException {
        int height = array.length;
        int width = array[0].length;

        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_BYTE_GRAY);

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                int intensity = array[y][x];
                if (intensity < 0 || intensity > 255) {
                    errorGeneral("A intensidade do pixel deve estar entre 0 e 255.");
                }
                int rgb = (intensity << 16) | (intensity << 8) | intensity;
                image.setRGB(x, y, rgb);
            }
        }

        File outputFile = new File(outputFilePath);
        ImageIO.write(image, "png", outputFile);
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
            file = new File(outputFolderPath + "/" + pngFileName.replace(".jpg", "(" + counter + ").jpg"));
            counter++;
        }

        try {
            writeArrayAsImage(normalizedImage, outputPath);
            if (printOrNot == 1) {
                System.out.println("A imagem mais proxima foi salva com sucesso: " + outputPath);
            }
        } catch (IOException e) {
            System.err.println("Erro ao salvar a imagem: " + e.getMessage());
        }
    }

    public static void saveMatrixToFile(double[][] matrix, String inputCsvPath, String outputFolderPath, int printOrNot) {
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

    public static String receiveImageLocation(String[] args) {
        String imageFolderLocationArgs;
        if (args == null) {
            String imageFolderLocation = scanner.next();
            if (!checkImagesFolderLocation(imageFolderLocation)) {
                System.out.println("Erro: Localização da base de imagens inválida");
                System.out.println("Tentar novamente ? (S/N)");
                String answer = scanner.next().toUpperCase();
                if (answer.equals("S")) {
                    imageFolderLocation = verifyCsvLocation();
                } else {
                    System.out.println("Saindo da aplicação, ainda pode desistir mas retornará ao menu inicial.");
                    quitApplication();
                }
            }
            return imageFolderLocation;
        } else {
            imageFolderLocationArgs = args[5];
            if (!checkImagesFolderLocation(imageFolderLocationArgs)) {
                errorGeneral("Erro: Localização inválida csv");
            }
            return imageFolderLocationArgs;
        }
    }

    public static String receiveCsvLocation(String[] args) {
        String csvLocationArgs;
        if (args == null) {
            String csvLocation = scanner.next();
            if (!checkCsvLocation(csvLocation)) {
                System.out.println("Erro: Localização inválida CSV");
                System.out.println("Tentar novamente? (S/N)");
                String answer = scanner.next().toUpperCase();
                if (answer.equals("S")) {
                    csvLocation = verifyCsvLocation();
                } else {
                    System.out.println("Saindo da aplicação, ainda pode desistir mas retornará ao menu inicial.");
                    quitApplication();
                }
            }

            csvLocation = verifySymmetricMatrix(csvLocation);

            return csvLocation;
        } else {
            csvLocationArgs = args[5];
            if (!checkCsvLocation(csvLocationArgs)) {
                errorGeneral("Erro: Localização inválida csv");
            }

            csvLocationArgs = verifySymmetricMatrix(csvLocationArgs);

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

    public static double[][] getMatrixFromCsv(String csvLocation) {
        int[] dimensions = getDimensions();
        int rows = dimensions[0];
        int cols = dimensions[1];

        double[][] matrix = new double[rows][cols];
        populateMatrix(matrix, csvLocation);

        return matrix;
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
    //* ------------------ Fim dos metodos de entrada e saída ------------------


    //* ------------------ Verificações ------------------
    public static boolean checkIfIsSymmetric(double[][] matrix) {
        int a = matrix.length;
        for (int i = 0; i < a; i++) {
            for (int j = i + 1; j < a; j++) {
                if (matrix[i][j] != matrix[j][i]) {
                    return false;
                }
            }
        }
        return true;
    }
    public static boolean checkCorrectParametersStructure(String[] parameters) {
        if (parameters.length == 8) {
            return parameters[0].equals("-f") && parameters[2].equals("-k") && parameters[4].equals("-i") && parameters[6].equals("-j");
        }
        return false;
    }

    public static boolean checkFunctionOptions(int function) {
        return function >= 1 && function <= 4;
    }

    public static boolean checkSizeBoundaries(int rows, int cols) {
        return rows > MAX_SIZE_ROWS || cols > MAX_SIZE_COLS || rows < MIN_SIZE_ROWS || cols < MIN_SIZE_COLS;
    }

    public static boolean checkCsvLocation(String csvLocation) {
        File csv = new File(csvLocation);
        if (csvLocation.isEmpty()) {
            return false;
        } else if (!csvLocation.contains(".csv")) {
            return false;
        } else return csv.exists();
    }

    public static boolean checkImagesFolderLocation(String imageFolderLocation) {
        File imageDirectory = new File(imageFolderLocation);
        if (imageFolderLocation.isEmpty()) {
            return false;
        }
        return imageDirectory.exists();
    }
    public static String verifySymmetricMatrix(String csvLocation) {

        double[][] matrix = readCSVToMatrix(csvLocation);

        while (!checkIfIsSymmetric(matrix)) {
            System.out.println("A matriz não é simétrica.");
            System.out.println("Tentar novamente? (S/N)");
            String answer = scanner.next().toUpperCase();
            if (answer.equals("S")) {
                csvLocation = verifyCsvLocation();
                matrix = readCSVToMatrix(csvLocation);
            } else {
                System.out.println("Saindo da aplicação, ainda pode desistir mas retornará ao menu inicial.");
                quitApplication();
            }
        }


        return csvLocation;
    }

    public static String verifyCsvLocation() {
        String csvLocation;
        uiCsvLocationParameterMenu();
        csvLocation = receiveCsvLocation(null);
        return csvLocation;
    }
    public static String verifyImageFolderLocation(int function) {
        String imageFolderLocation;
        do {
            uiImageLocationParameterMenu();
            imageFolderLocation = receiveImageLocation(null);
        } while (!checkImagesFolderLocation(imageFolderLocation));
        return imageFolderLocation;
    }

    public static int verifyFunction() {
        int function;
        do {
            uiFunctionParameterMenu();
            function = receiveFunction(null);
        } while (!checkFunctionOptions(function));
        return function;
    }

    public static int verifyVectorNumbers() {
        int vectorNumbers;
        do {
            uiVectorNumbersParameterMenu();
            vectorNumbers = receiveNumberVectors(null);
        } while (vectorNumbers <= MIN_QUANTITY_VECTORS);
        return vectorNumbers;
    }

    public static void checkExistanceFileDirectory(String csvLocation) {
        try {
            scannerCsv = new Scanner(new File(csvLocation));
        } catch (FileNotFoundException e) {
            errorGeneral("Erro ao abrir os arquivos: " + e.getMessage());
        }
    }
    //* ------------------ Fim verificações ------------------


    //* ------------------ Menus de opções ------------------
    public static void uiFunctionParameterMenu() {
        System.out.println(" ____________________________________________________");
        System.out.println("|            Que função deseja realizar?             |");
        System.out.println("|                                                    |");
        System.out.println("| 1 - Decomposição Própria de uma Matriz Simétrica.  |");
        System.out.println("| 2 - Reconstrução de Imagens usando Eigenfaces.     |");
        System.out.println("| 3 - Identificação de imagem mais próxima.          |");
        System.out.println("| 4 - Deseja sair da aplicação ?                     |");
        System.out.println(" ----------------------------------------------------");
        System.out.print("Opção: ");
    }

    public static void uiVectorNumbersParameterMenu() {
        System.out.println("----- Quantos vetores próprios deseja utilizar? -----");
        System.out.print("Quantidade: ");
    }

    public static void uiCsvLocationParameterMenu() {
        System.out.println("-- Qual a localização do csv que deseja utilizar?  --");
        System.out.print("Localização: ");
    }

    public static void uiImageLocationParameterMenu() {
        System.out.println("------- Qual a localização da base de imagens? ------");
        System.out.print("Localização: ");
    }

    public static void uiQuitParameterMenu() {
        System.out.println("|-- Tem certeza que deseja sair da aplicação? (S/N) --|");
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
                errorGeneral("Erro: Opção inválida");
            }
            return functionArgs;
        }
    }

    public static int receiveNumberVectors(String[] args) {
        int vectorNumbersArgs;
        if (args == null) {
            vectorNumbersArgs = scanner.nextInt();
            if (vectorNumbersArgs <= MIN_QUANTITY_VECTORS) {
                System.out.println("Erro: O número de vetores deve ser maior que 0.");
            }
        } else {
            vectorNumbersArgs = Integer.parseInt(args[3]);
            if (vectorNumbersArgs <= MIN_QUANTITY_VECTORS) {
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
                    System.out.println("Retornando ao menu de opções.");
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
    //* ----------------- Fim operaçoes básicas com matrizes ------------------


    //* -------------------- Printar Matrizes -----------------------
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

    public static void printLine(int length, String pattern) {
        for (int i = 0; i < length; i++) {
            System.out.print(pattern);
        }
        System.out.println();
    }
    //* ----------------- Fim printar matrizes -----------------------


    //* ----------------- Printar Funcionalidades -----------------
    public static void printFunction1(int numbersEigenfaces, double[][] newEigenValuesK, double[][] newEigenVectorsK, double maximumAbsolutError) {
        System.out.println("A quantidade selecionada para a variável K foi: " + numbersEigenfaces);
        printMatrix(newEigenValuesK, "Matriz Valores Próprios K");
        printMatrix(newEigenVectorsK, "Matriz Vetores Próprios K:");
        System.out.printf("Erro Absoluto Médio: %.3f\n", maximumAbsolutError);
    }

    public static void printHeaderFunction(String functionName) {
        System.out.println();
        printLine(1, "-------------------------------------------------------");
        System.out.println(functionName);
        printLine(1, "-------------------------------------------------------");
        System.out.println();
    }
    //* ----------------- Fim printar funcionalidades -----------------------


    //! ------------------ Error Messages ------------------
    public static void errorGeneral(String error) {
        System.out.println(error);
        System.exit(1);
    }
    //! ------------------ Fim error messages --------------
}