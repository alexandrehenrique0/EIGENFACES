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

//! PRECISA SER ESTUDADO MAS JÀ ESTÁ FUNCIONANDO PARA VARIAS IMAGENS
public class GabrielTesteFunc2 {
    public static void main(String[] args) {
        String csvPath = "Input/Funcao2-3/csv";

        // Lista os nomes dos arquivos CSV
        String[] CSVFileNames = getCSVFileNames(csvPath);

        // Carrega todas as matrizes
        double[][][] matrixCSVDouble3D = get_Matrices_From_CSV_Folder(csvPath);

        // Consolidar todas as imagens em uma matriz 2D (linearizada)
        double[][] linearizedImages = new double[matrixCSVDouble3D[0].length * matrixCSVDouble3D[0][0].length][matrixCSVDouble3D.length];
        for (int img = 0; img < matrixCSVDouble3D.length; img++) {
            linearizedImages = flattenImagesMatrix(matrixCSVDouble3D[img]);
        }

        // Calcula o vetor médio global
        double[] meanVector = calculateMeanVector(linearizedImages);

        // Centraliza todas as imagens (Matriz Φ)
        double[][] phi = centralizeImages(linearizedImages, meanVector);

        // Calcula valores e vetores próprios com base na matriz centralizada
        int numOfColumns = matrixCSVDouble3D[0][0].length; // Número de colunas de uma matriz original
        double[][] covariance = covariances(phi, numOfColumns);
        double[][] eigenVectors = eigenVectors(covariance, phi);
        double[][] normalizedVectors = normalizeVectors(eigenVectors);

        // Itera sobre cada imagem para calcular pesos e reconstruir
        for (int img = 0; img < linearizedImages[0].length; img++) {
            // Projeta a imagem atual nos eigenfaces
            double[] weights = calculateWeights(phi[img], normalizedVectors);

            // Reconstrói a imagem
            double[] reconstructedImage = reconstructImage(meanVector, normalizedVectors, weights, numOfColumns);

            // Converte para matriz 2D
            double[][] reconstructedImageMatrix = array1DToMatrix(reconstructedImage, matrixCSVDouble3D[img]);

            // Ajusta a precisão e salva a imagem
            adjustPrecision(reconstructedImageMatrix);
            saveImage(reconstructedImageMatrix, CSVFileNames[img], "Output/Func2/ImagensReconstruidas");
            saveMatrixToFile(reconstructedImageMatrix, CSVFileNames[img], "Output/Func2/Eigenfaces");
        }


        //! Print para verificar cada etapa
//        print_Matrix(originalMatrix, "originalMatrix");
//        System.out.println("Número de colunas: " + numOfColumns);
//        System.out.println("Número de eigenfaces: " + eigenFaces);
//        print_Matrix(averageColumn, "averageColumn");
//        print_Matrix(deviantMatrix, "deviantMatrix");
//        print_Matrix(deviantMatrixTranspose, "deviantMatrixTranspose");
//        print_Matrix(covariance, "covariance");
//        print_Matrix(eigenValuesATxA, "eigenValuesATxA");
//        print_Matrix(eigenVectorsATxA, "eigenVectorsATxA");
//        print_Matrix(eigenVectorsAxAT, "eigenVectorsAxAT");
//        print_Matrix(eigenValuesAxAT, "eigenValuesAxAT");
//        print_Matrix(normalizedVectorsATxA, "normalizedVectorsATxA");
//        System.out.println("matrixLinearized lenght: " + matrixLinearized.length);
//        System.out.println("matrixLinearized: " + Arrays.toString(matrixLinearized));
//        System.out.println("linearizedImages lenght: " + linearizedImages.length);
//        print_Matrix(linearizedImages, "linearizedImages");
//        System.out.println("meanVector lenght: " + meanVector.length);
//        System.out.println("meanVector: " + Arrays.toString(meanVector));
//        System.out.println("phi lenght: " + phi.length);
//        System.out.println("phi: " + Arrays.toString(phi));
//        System.out.println("weights lenght: " + weights.length);
//        System.out.println("weights: " + Arrays.toString(weights));
//        System.out.println("reconstructedImage lenght: " + reconstructedImage.length);
//        System.out.println("reconstructedImage: " + Arrays.toString(reconstructedImage));
//        print_Matrix(reconstructedImageMatrix, "reconstructedImageMatrix");

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

        int index = 0;
        for (File csvFile : csvFiles) {
            matrices[index] = readCSVToArray(csvFile.getPath());
            index++;
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
                String rowString = String.join(",", Arrays.stream(row)
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

    private static double[][] flattenImagesMatrix(double[][] matrixOfImages) {
        int rows = matrixOfImages.length;
        int columns = matrixOfImages[0].length;

        double[][] linearizedImages = new double[rows * columns][1];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                linearizedImages[i * columns + j][0] = matrixOfImages[i][j];
            }
        }
        return linearizedImages;
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
        String pngFileName = csvFileName.replace(".csv", ".png");
        String outputPath = outputFolderPath + "/" + pngFileName;

        try {
            writeArrayAsImage(normalizedImage, outputPath);
            System.out.println("Image saved successfully: " + outputPath);
        } catch (IOException e) {
            System.err.println("Failed to save the image: " + e.getMessage());
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


    //! se para calcular os pesos pega phi e eigenfaces, revisar pois esta recebendo uma matriz
    //! de vetores normalizados e não a quantidade de eigenfaces (que nesse caso é a mesma qtd de colunas)
    public static double[] calculateWeights(double[] phi, double[][] eigenfaces) {
        if (Math.sqrt(phi.length) != eigenfaces.length) {
            throw new IllegalArgumentException("O comprimento de 'phi' deve ser igual ao número de linhas em 'eigenfaces'.");
        }
        double[] weights = new double[eigenfaces[0].length];
        for (int j = 0; j < eigenfaces[0].length; j++) {
            for (int i = 0; i < eigenfaces[0].length; i++) {
                weights[j] += phi[i] * eigenfaces[i][j];
            }
        }
        return weights;
    }

    //! revisar essa matriz que recebe chamada eigenFaces
    public static double[] reconstructImage(double[] meanVector, double[][] eigenfaces, double[] weights, int quantityEigenfaces) {
        quantityEigenfaces = Math.min(quantityEigenfaces, eigenfaces[0].length);

        double[] reconstructed = new double[meanVector.length];
        for (int i = 0; i < meanVector.length; i++) {
            reconstructed[i] = meanVector[i];
        }

        for (int j = 0; j < quantityEigenfaces; j++) {
            for (int i = 0; i < eigenfaces[0].length; i++) {
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

    private static void print_Line(int length, String pattern) {
        for (int i = 0; i < length; i++) {
            System.out.print(pattern);
        }
        System.out.println();
    }

    private static void print_Matrix(double[][] matrixToPrint, String matrixName) {
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

    //* calcula a media de cada coluna da matriz
    public static double[][] averageColumn(double[][] matrix) {
        int rows = matrix.length;
        int columns = matrix[0].length;
        double[][] averageColumn = new double[rows][1];
        for (int i = 0; i < rows; i++) {
            double sum = 0;
            for (int j = 0; j < columns; j++) {
                sum = sum + matrix[i][j];
            }
            averageColumn[i][0] = sum / columns;
        }
        return averageColumn;
    }

    public static double[][] deviantMatrix(double[][] matrix, double[][] averageColumn) {
        int columns = matrix[0].length;
        int rows = matrix.length;
        double[][] deviant = new double[rows][columns];

        for (int j = 0; j < columns; j++) {
            for (int i = 0; i < rows; i++) {
                deviant[i][j] = deviantColumn(matrix[i][j], averageColumn[i][0]); //! averageColumn[i][0] é a media de cada linha
            }

        }
        return deviant;
    }

    public static double deviantColumn(double matrixValue, double averageColumnValue) {

        double deviant = matrixValue - averageColumnValue; // Calcula o desvio de cada elemento da matriz em relação à média da linha correspondente

        return deviant;
    }

    public static double[][] covariances(double[][] matrixA, int quantityOfColumns) {
        double[][] matrixATransposed = transposeMatrix(matrixA);
        double[][] matrixAmultplyByAT = multiplyMatrix(matrixA, matrixATransposed);
        return multiplyMatrixByScalar(matrixAmultplyByAT, 1.0 / quantityOfColumns);
    }

    public static double[][] eigenVectors(double[][] matrixOne, double[][] matrixTwo) {
        double[][] matrixThree = multiplyMatrix(matrixOne, matrixTwo);
        EigenDecomposition decomposedMatrixThree = decomposeMatrix(matrixThree);
        RealMatrix eigenVectorsMatrix = decomposedMatrixThree.getV();
        return eigenVectorsMatrix.getData();
    }

    public static double[][] eigenValues(double[][] matrixOne, double[][] matrixTwo) {
        double[][] matrixThree = multiplyMatrix(matrixOne, matrixTwo);
        EigenDecomposition decomposedMatrixThree = decomposeMatrix(matrixThree);
        RealMatrix eigenVectorsMatrix = decomposedMatrixThree.getD();
        return eigenVectorsMatrix.getData();
    }

    //ignora valores negativos muito baixos, para nao obter -0.0
    private static void adjustPrecision(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (Math.abs(matrix[i][j]) < 0.001) {
                    matrix[i][j] = 0.0;
                }
            }
        }
    }

    public static double[][] normalizeVectors(double[][] eigenVectorsATxA) {
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
        double[][] resultMatrix = new double[matrixLeft.length][matrixRight[0].length];
        for (int i = 0; i < matrixLeft.length; i++) {
            for (int j = 0; j < matrixRight[0].length; j++) {
                for (int k = 0; k < matrixRight.length; k++) {
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

//public static void main(String[] args) {
//        String csvPath = "Input/Funcao2-3/csv";
//
//        String[] CSVFileNames = getCSVFileNames(csvPath);
//
//        double[][][] matrixCSVDouble3D = get_Matrices_From_CSV_Folder(csvPath);
//
//        for (int i = 0; i < matrixCSVDouble3D.length; i++) {
//            double[][] originalMatrix = matrixCSVDouble3D[i];
//            String outputPath = CSVFileNames[i];
//
//            // Dimensão e número de eigenfaces
//            int numOfColumns = originalMatrix[0].length;
//            int eigenFaces = numOfColumns; // ajuste conforme necessário, mas se tiver colunas = eigenfaces imagem fica melhor
//
//            //* 1. Cálculo da média - revisto
//            double[][] averageColumn = averageColumn(originalMatrix);
//            adjustPrecision(averageColumn);
//
//            //* 2. Centralização (Matriz de desvios) - revisto
//            double[][] deviantMatrix = deviantMatrix(originalMatrix, averageColumn);
//            adjustPrecision(deviantMatrix);
//
//            double[][] deviantMatrixTranspose = transposeMatrix(deviantMatrix);
//            adjustPrecision(deviantMatrixTranspose);
//
//            //* 3. Matriz de Covariância - revisto
//            double[][] covariance = covariances(deviantMatrix, numOfColumns);
//            adjustPrecision(covariance);
//
//            //* 4. Valores e vetores próprios - revisto
//            double[][] eigenValuesATxA = eigenValues(deviantMatrixTranspose, deviantMatrix);
//            adjustPrecision(eigenValuesATxA);
//            double[][] eigenVectorsATxA = eigenVectors(deviantMatrixTranspose, deviantMatrix);
//            adjustPrecision(eigenVectorsATxA);
//
//            double[][] eigenValuesAxAT = eigenValues(deviantMatrix, deviantMatrixTranspose);
//            adjustPrecision(eigenValuesAxAT);
//            double[][] eigenVectorsAxAT = eigenVectors(deviantMatrix, deviantMatrixTranspose);
//            adjustPrecision(eigenVectorsAxAT);
//
//            double[][] normalizedVectorsATxA = normalizeVectors(eigenVectorsATxA);
//            adjustPrecision(normalizedVectorsATxA);
//
//            //* 5. Centralização da imagem para reconstrução - revisto
//            double[] matrixLinearized = matrixToArray1D(originalMatrix);
//
//            //! futuramente para várias imagens
//            double[][] linearizedImages = flattenImagesMatrix(originalMatrix);
//            //! precisará adaptar o código daqui pra baixo para suportar várias imagens
//
//            //* Calcula o vetor médio para o vetor linearizado - revisto
//            double[] meanVector = calculateMeanVector(linearizedImages);
//
//            //* Centraliza a imagem original linearizada - revisto
//            //! pois depende de meanVector que precisa ser revisto
//            double[] phi = centralizeImage(matrixLinearized, meanVector);
//
//            //* 6. Cálculo dos pesos (projeção nos eigenfaces) - revisto
//            double[] weights = calculateWeights(phi, normalizedVectorsATxA);
//
//            //* 7. Reconstrução da imagem com eigenfaces - revisto
//            double[] reconstructedImage = reconstructImage(meanVector, normalizedVectorsATxA, weights, eigenFaces);
//
//            //* 8. Conversão para matriz 2D
//            double[][] reconstructedImageMatrix = array1DToMatrix(reconstructedImage, originalMatrix);
//            adjustPrecision(reconstructedImageMatrix);
//
//            //* 9. Salvar a imagem reconstruída
//            saveImage(reconstructedImageMatrix, outputPath, "Output/Func2/ImagensReconstruidas");
//            saveMatrixToFile(reconstructedImageMatrix, outputPath, "Output/Func2/Eigenfaces");
//        }
