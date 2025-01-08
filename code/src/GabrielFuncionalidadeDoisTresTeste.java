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

//! classe para testar func 2 e 3
public class GabrielFuncionalidadeDoisTresTeste {
    public static void main(String[] args) {
        String csvPath = "Input/Funcao2-3/csv";

        String[] CSVFileNames = getCSVFileNames(csvPath);

        double[][] matrixCSVDouble = readCSVToArray("Input/Funcao2-3/csv/image_271.csv");

        double[][][] matrixCSVDouble3D = get_Matrices_From_CSV_Folder(csvPath);

        double[][] linearizedImages = new double[matrixCSVDouble3D[0].length * matrixCSVDouble3D[0].length][matrixCSVDouble3D.length];
        for (int img = 0; img < matrixCSVDouble3D.length; img++) {
            double[] matrixLinearized = matrixToArray1D(matrixCSVDouble3D[img]);
            for (int i = 0; i < matrixLinearized.length; i++) {
                linearizedImages[i][img] = matrixLinearized[i];
            }
        }

        double[] meanVector = calculateMeanVector(linearizedImages);

        double[][] phi = centralizeImages(linearizedImages, meanVector);

        double[][] phiT = transposeMatrix(phi);

        double[][] ATxA = multiplyMatrix(phiT, phi);

        double[][] eigenVectors = eigenVectors(ATxA);

        double[][] eigenVectorsExpand = multiplyMatrix(phi, eigenVectors);

        double[][] eigenfaces = normalize(eigenVectorsExpand);

        double[][] weightsMatrix = new double[phi.length][phi[0].length];

        for (int img = 0; img < linearizedImages[0].length; img++) {
            double[] weights = calculateWeightsOne(getColumn(phi, img), eigenfaces);

            for (int i = 0; i < weights.length; i++) {
                weightsMatrix[i][img] = weights[i];
            }

            System.out.println("pesos");
            System.out.println(Arrays.toString(weights));

            double[] reconstructedImage = reconstructImage(meanVector, eigenfaces, weights);

            double[][] reconstructedImageMatrix = array1DToMatrix(reconstructedImage, matrixCSVDouble3D[img]);

            saveImage(reconstructedImageMatrix, CSVFileNames[img], "Output/Func2/ImagensReconstruidas");
            saveMatrixToFile(reconstructedImageMatrix, CSVFileNames[img], "Output/Func2/Eigenfaces");
        }

//        double[] phiVector = centralizeVector(matrixToArray1D(matrixCSVDouble), meanVector);
//        double[] weightsVetorPrincipal = calculateWeights(phiVector, transposeMatrix(normalizedEigenfaces));
//
//        double[] matrizEuclidiana = calculate_Euclidian_Distance(weightsVetorPrincipal, weightsMatrix);
//        int posicaoMaisProxima = check_Closer_Vetor(matrizEuclidiana);

//        System.out.println("A imagem mais próxima é: " + CSVFileNames[posicaoMaisProxima]);
    }

    private static void adjustPrecision(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (Math.abs(matrix[i][j]) < 1e-10) { // Limite para evitar valores muito pequenos
                    matrix[i][j] = 0.0;
                }
            }
        }
    }

    public static void printMatrix(double[][] matrixToPrint, String matrixName) {
        System.out.println("Matriz: " + matrixName + " ↓");
        printLine(matrixToPrint[0].length, "____________");

        for (int j = 0; j < matrixToPrint.length; j++) {
            double[] row = matrixToPrint[j];
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
                        .mapToObj(val -> String.format("%.0f", val))
                        .toArray(String[]::new));
                writer.println(rowString);
            }
            System.out.println("Arquivo CSV criado com sucesso: " + file.getName());
        } catch (IOException e) {
            System.err.println("Erro ao salvar a matriz no arquivo: " + e.getMessage());
            e.printStackTrace();
        }
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
            file = new File(outputFolderPath + "/" + pngFileName.replace(".jpg", "(" + counter + ").jpg"));
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
                if (array[i * columns + j] < 0) {
                    matrix[i][j] = 0;
                } else if (array[i * columns + j] > 255) {
                    matrix[i][j] = 255;
                } else {
                    matrix[i][j] = array[i * columns + j];
                }
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
                if (intensity < 0) {
                    intensity = 0;
                }
                if (intensity > 255){
                    intensity = 255;
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
        double[][] matrixATmultiplyByA = multiplyMatrix(matrixATransposed, matrixA);
        return multiplyMatrixByScalar(matrixATmultiplyByA, 1.0 / quantityOfImages);
    }

    public static double[][] eigenVectors(double[][] matrix) {
        EigenDecomposition decomposedMatrixThree = decomposeMatrix(matrix);
        RealMatrix eigenVectorsMatrix = decomposedMatrixThree.getV();
        return eigenVectorsMatrix.getData();
    }

    public static double[][] normalize(double[][] eigenVectorsATxA) {
        double[][] eigenVectorsATxACopy = new double[eigenVectorsATxA.length][eigenVectorsATxA[0].length];

        for (int i = 0; i < eigenVectorsATxA.length; i++) {
            for (int j = 0; j < eigenVectorsATxA[0].length; j++) {
                eigenVectorsATxACopy[i][j] = eigenVectorsATxA[i][j];
            }
        }

        for (int i = 0; i < eigenVectorsATxACopy[0].length; i++) {
            double norm = 0;

            for (int j = 0; j < eigenVectorsATxACopy.length; j++) {
                norm += eigenVectorsATxACopy[j][i] * eigenVectorsATxACopy[j][i];
            }
            norm = Math.sqrt(norm);

            for (int j = 0; j < eigenVectorsATxA.length; j++) {
                if (norm == 0) {
                    eigenVectorsATxACopy[j][i] = 0;
                } else {
                    eigenVectorsATxACopy[j][i] /= norm;
                }
            }
        }
        return eigenVectorsATxACopy;
    }

    public static EigenDecomposition decomposeMatrix(double[][] matrix) {
        Array2DRowRealMatrix matrixToDecompose = new Array2DRowRealMatrix(matrix);
        return new EigenDecomposition(matrixToDecompose);
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

    public static double[][] multiplyMatrix(double[][] matrixLeft, double[][] matrixRight) {
        if (matrixLeft[0].length != matrixRight.length) {
            throw new IllegalArgumentException("Matrizes incompatíveis para multiplicação: colunas de matrixLeft != linhas de matrixRight");
        }

        double[][] resultMatrix = new double[matrixLeft.length][matrixRight[0].length];
        for (int i = 0; i < matrixLeft.length; i++) {
            for (int j = 0; j < matrixRight[0].length; j++) {
                resultMatrix[i][j] = 0;
                for (int k = 0; k < matrixRight.length; k++) {
                    double partialProduct = matrixLeft[i][k] * matrixRight[k][j];
                    if (Double.isNaN(partialProduct) || Double.isInfinite(partialProduct)) {
                        System.out.println("Erro: Valor inválido em matrixLeft[" + i + "][" + k + "] * matrixRight[" + k + "][" + j + "]");
                        System.out.println("matrixLeft[" + i + "][" + k + "] = " + matrixLeft[i][k]);
                        System.out.println("matrixRight[" + k + "][" + j + "] = " + matrixRight[k][j]);
                        System.exit(1);
                    }
                    resultMatrix[i][j] += matrixLeft[i][k] * matrixRight[k][j];
                }
            }
        }
        return resultMatrix;
    }

    public static double[][] multiplyMatrixByScalar(double[][] matrix, double scalar) {
        double[][] resultMatrix = new double[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                resultMatrix[i][j] = matrix[i][j] * scalar;
            }
        }
        return resultMatrix;
    }
}
