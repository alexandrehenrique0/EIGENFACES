import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class GabrielTesteFunc2 {
    public static void main(String[] args) {
        String directoryPath = "Input/TesteFuncao2-3/csv";
        File dir = new File(directoryPath);
        //! Adicionado log para verificar se o diretório existe e contém arquivos
        if (!dir.exists() || !dir.isDirectory()) {
            System.out.println("Erro: O diretório especificado não existe ou não é um diretório válido.");
            System.exit(1);
        }

        File[] csvFiles = dir.listFiles((d, name) -> name.endsWith(".csv"));
        //! Validando se os arquivos CSV foram encontrados
        if (csvFiles == null || csvFiles.length == 0) {
            System.out.println("Nenhuma imagem encontrada no diretório.");
            System.exit(1);
        }

        //! Ordenar arquivos para consistência, caso necessário
        Arrays.sort(csvFiles);

        String csvPath = csvFiles[0].getAbsolutePath();
        System.out.println("Primeiro arquivo CSV selecionado: " + csvPath);

        double[][] firstImage = readCSVToArray(csvFiles[0].getAbsolutePath());
        System.out.println("Dimensões da primeira imagem: " + firstImage.length + "x" + firstImage[0].length);

        int M = firstImage.length;
        int N = firstImage[0].length;
        int numPixels = M * N;

        int numImagens = csvFiles.length;
        System.out.println("Número de imagens encontradas: " + numImagens);

        double[][] originalArray = new double[numPixels][numImagens];

        //! Log ao carregar imagens
        System.out.println("Carregando imagens na matriz total...");
        loadImagesToMatrix(csvFiles, originalArray, firstImage);
        System.out.println("Todas as imagens carregadas na matriz.");

        //! Calculando vetor médio
        System.out.println("Calculando vetor médio...");
        double[] meanVectors = meanVectors(numPixels, numImagens, originalArray);
        System.out.println("Vetor médio calculado com sucesso.");

        double[][] original = readCSVToArray(csvPath);
        System.out.println("Matriz original carregada do primeiro CSV.");

        //* Para reconstrução é interessante que os Eigenvectors sejam iguais à quantidade de colunas da matriz original
        //! Ajuste para garantir que `numEigenfaces` não ultrapasse o limite das colunas
        int numEigenfaces = Math.min(originalArray[0].length, original[0].length);
        System.out.println("Número de Eigenfaces definido: " + numEigenfaces);

        //* 1. Cálculo da média - revisto
        double[] averageColumn = averageColoumn(originalArray);
        System.out.println("Média das colunas calculada.");

        //* 2. Centralização (Matriz de desvios) - revisto
        double[][] deviationMatrixA = deviationMatrix(originalArray, averageColumn);
        System.out.println("Dimensões de deviationMatrixA: " + deviationMatrixA.length + "x" + deviationMatrixA[0].length);
        adjustPrecision(deviationMatrixA);

        double[][] deviationMatrixATranspose = transposeMatriz(deviationMatrixA);
        adjustPrecision(deviationMatrixATranspose);
        System.out.println("Matriz de desvios e transposta calculadas e ajustadas.");

        //* 3. Matriz de Covariância - revisto
        double[][] covariance = covariances(deviationMatrixA);
        adjustPrecision(covariance);
        System.out.println("Matriz de covariância calculada.");

        //* Buscar matriz AxAT e ATxA e decompor - revisto
        double[][] ATxA = multiplyMatrix(deviationMatrixATranspose, deviationMatrixA);
        adjustPrecision(ATxA);
        EigenDecomposition eigenDecompositionATxA = decomposeMatrix(ATxA);

        double[][] AxAT = multiplyMatrix(deviationMatrixA, deviationMatrixATranspose);
        adjustPrecision(AxAT);
        EigenDecomposition eigenDecompositionAxAT = decomposeMatrix(AxAT);
        System.out.println("Decomposição de matrizes (ATxA e AxAT) concluída.");

        //* 4. Valores e vetores próprios - revisto
        double[][] eigenValuesATxA = getEigenValues(eigenDecompositionATxA);
        adjustPrecision(eigenValuesATxA);
        double[][] eigenVectorsATxA = getEigenVectors(eigenDecompositionATxA);
        adjustPrecision(eigenVectorsATxA);

        double[][] eigenVectorsAxAT = getEigenVectors(eigenDecompositionAxAT);
        adjustPrecision(eigenVectorsAxAT);
        System.out.println("Dimensões de eigenVectorsAxAT: " + eigenVectorsAxAT.length + "x" + eigenVectorsAxAT[0].length);

        double[][] eigenValuesAxAT = getEigenValues(eigenDecompositionAxAT);
        adjustPrecision(eigenValuesAxAT);
        System.out.println("Valores e vetores próprios calculados e ajustados.");



        double[][] normalizedEigenVectors = normalizeVectors(eigenVectorsAxAT); // Normalizar os vetores próprios compatíveis
        adjustPrecision(normalizedEigenVectors);
        System.out.println("Dimensões de normalizedEigenVectors: " + normalizedEigenVectors.length + "x" + normalizedEigenVectors[0].length);

        // Reduzir normalizedEigenVectors para os 40 principais Eigenfaces
        double[][] reducedEigenVectors = new double[normalizedEigenVectors.length][numEigenfaces];
        for (int i = 0; i < normalizedEigenVectors.length; i++) {
            for (int j = 0; j < numEigenfaces; j++) {
                reducedEigenVectors[i][j] = normalizedEigenVectors[i][j];
            }
        }
        normalizedEigenVectors = reducedEigenVectors;
        System.out.println("Dimensões de normalizedEigenVectors após redução: " + normalizedEigenVectors.length + "x" + normalizedEigenVectors[0].length);


        //* 5. Centralização da imagem para reconstrução - revisto
        double[] linearizedOriginalMatrix = matrixToArray(original);
        System.out.println("Tamanho de linearizedOriginalMatrix: " + linearizedOriginalMatrix.length);
        System.out.println("Tamanho de meanVectors: " + meanVectors.length);

        //* Centraliza a imagem original linearizada - revisto
        if (linearizedOriginalMatrix.length != meanVectors.length) {
            System.out.println("Erro: Os tamanhos dos arrays 'linearizedOriginalMatrix' e 'meanVector' não correspondem.");
            System.exit(1);
        }
        double[] phi = centralizeImagePhi(linearizedOriginalMatrix, meanVectors);
        System.out.println("Imagem original centralizada (phi) calculada.");
        //* Validando se as dimensões de 'normalizedEigenVectors' correspondem ao tamanho de 'phi'
        if (normalizedEigenVectors.length != phi.length) {
            System.out.println("Erro: O número de linhas em 'normalizedEigenVectors' não corresponde ao tamanho de 'phi'.");
            System.out.println("Dimensões de normalizedEigenVectors: " + normalizedEigenVectors.length + "x" + normalizedEigenVectors[0].length);
            System.out.println("Tamanho de phi: " + phi.length);
            System.exit(1);
        }

        //* 6. Cálculo dos pesos (projeção nos eigenfaces) - revisto
        double[] weights = calculateWeights(phi, normalizedEigenVectors);
        System.out.println("Tamanho de weights: " + weights.length);

        //* 7. Reconstrução da imagem com eigenfaces (1D) - revisto
        double[] reconstructedImageArray = reconstructImage(averageColumn, normalizedEigenVectors, weights, numEigenfaces);
        System.out.println("Imagem reconstruída (1D).");


        System.out.println("Dimensões de normalizedEigenVectors: " + normalizedEigenVectors.length + "x" + normalizedEigenVectors[0].length);
        System.out.println("Tamanho de weights: " + weights.length);
        System.out.println("Número de Eigenfaces usadas: " + numEigenfaces);

        //* 8. Conversão para matriz 2D - revisto
        // Reconstruir a matriz com as dimensões da imagem original (64x64)
        double[][] reconstructedImageMatrix = arrayToMatrix(reconstructedImageArray, new double[64][64]);
        adjustPrecision(reconstructedImageMatrix);
        System.out.println("Imagem reconstruída convertida para matriz 2D com dimensões 64x64.");

        //* 9. Salvar a imagem reconstruída - revisto
        saveImage(reconstructedImageMatrix, csvPath, "Input/TesteFuncao2-3/OutputImagesFunc2-3");
        System.out.println("Imagem reconstruída salva com sucesso.");
    }


    private static double[][] loadImagesToMatrix(File[] csvFiles, double[][] matrixTotalImages, double[][] firstImage) {
        for (int idx = 0; idx < csvFiles.length; idx++) {
            double[][] img = readCSVToArray(csvFiles[idx].getAbsolutePath());
            // Verifica se todas têm o mesmo tamanho
            if (img.length != firstImage.length || img[0].length != firstImage[0].length) {
                System.out.println("As imagens não têm as mesmas dimensões.");
                System.exit(1);
            }
            double[] linearized = matrixToArray(img);
            for (int i = 0; i < firstImage[0].length * firstImage.length; i++) {
                matrixTotalImages[i][idx] = linearized[i];
            }
        }
        return matrixTotalImages;
    }

    public static double[] meanVectors (int numPixels, int numImagens, double[][] A) {
        double[] meanVector = new double[numPixels];
        for (int i = 0; i < numPixels; i++) {
            double sum = 0;
            for (int j = 0; j < numImagens; j++) {
                sum += A[i][j];
            }
            meanVector[i] = sum / numImagens;
        }
        return meanVector;
    }

    public static double[] matrixToArray(double[][] matrix) {
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

    public static void saveImage(double[][] imageArray, String inputCsvPath, String outputFolderPath) {
        int height = imageArray.length;
        int width = imageArray[0].length;

        // Normalize the image to the range 0-255
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;

        // Find min and max values
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

        // Create output file name
        String csvFileName = new File(inputCsvPath).getName();
        String pngFileName = csvFileName.replace(".csv", ".png");
        String outputPath = outputFolderPath + "/" + pngFileName;

        // Write the image
        try {
            writeArrayAsImage(normalizedImage, outputPath);
            System.out.println("Image saved successfully: " + outputPath);
        } catch (IOException e) {
            System.err.println("Failed to save the image: " + e.getMessage());
        }
    }

    public static double[] centralizeImagePhi (double[] linearizedOriginalMatrix, double[] meanVector) {
        if (linearizedOriginalMatrix.length != meanVector.length){
            System.out.println("Os tamanhos dos arrays 'linearizedOriginalMatrix' e 'meanVector' não correspondem.");
            System.exit(1); //! termina o programa, pois a função não pode continuar
        }
        double[] phi = new double[linearizedOriginalMatrix.length];

        for (int i = 0; i < linearizedOriginalMatrix.length; i++) {
            phi[i] = linearizedOriginalMatrix[i] - meanVector[i];
        }
        return phi;
    }

    public static double[] calculateWeights(double[] phi, double[][] normalizedEigenVectors) {
        int numEigenfaces = normalizedEigenVectors[0].length; // Número de Eigenfaces disponíveis

        if (phi.length != normalizedEigenVectors.length) {
            System.out.println("Erro: O comprimento de 'phi' deve ser igual ao número de linhas em 'eigenfaces'.");
            System.out.println("Tamanho de phi: " + phi.length);
            System.out.println("Número de linhas em 'normalizedEigenVectors': " + normalizedEigenVectors.length);
            System.exit(1); //! Termina o programa
        }

        // Inicializa o vetor de pesos com base no número de Eigenfaces a serem usados
        double[] weights = new double[numEigenfaces];

        // Calcula os pesos projetando 'phi' nos Eigenfaces
        for (int j = 0; j < numEigenfaces; j++) { // Itera sobre as colunas dos Eigenfaces
            for (int i = 0; i < normalizedEigenVectors.length; i++) { // Itera sobre as linhas (pixels)
                weights[j] += phi[i] * normalizedEigenVectors[i][j];
            }
        }

        return weights;
    }



    public static double[] reconstructImage(double[] meanVector, double[][] normalizedEigenVectors, double[] weights, int numEigenfaces) {
        // Limita o número de Eigenfaces ao menor valor entre o fornecido e o número disponível
        numEigenfaces = Math.min(numEigenfaces, normalizedEigenVectors[0].length);

        // Inicializa o vetor reconstruído com o vetor médio
        double[] reconstructedImageArray = new double[meanVector.length];
        for (int i = 0; i < meanVector.length; i++) {
            reconstructedImageArray[i] = meanVector[i];
        }

        // Adiciona a contribuição de cada Eigenface
        for (int j = 0; j < numEigenfaces; j++) { // Itera sobre os Eigenfaces utilizados (colunas)
            for (int i = 0; i < meanVector.length; i++) { // Itera sobre os pixels (linhas)
                reconstructedImageArray[i] += weights[j] * normalizedEigenVectors[i][j];
            }
        }
        return reconstructedImageArray;
    }




    public static double[][] arrayToMatrix(double[] array, double[][] originalMatrix) {
        int rows = originalMatrix.length;    // Linhas da matriz
        int columns = originalMatrix[0].length; // Colunas da matriz

        if (array.length != rows * columns) {
            System.out.println("Erro: O tamanho do vetor não corresponde às dimensões da matriz.");
            System.out.println("Tamanho do vetor: " + array.length);
            System.out.println("Dimensões esperadas da matriz: " + rows + "x" + columns);
            System.exit(1);
        }

        double[][] matrix = new double[rows][columns];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                matrix[i][j] = array[i * columns + j];
            }
        }
        return matrix;
    }



    //! metodo fornecido no moodle
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

    //! usar metodo já fornecido no main, acho que rafael já fez isso
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

    public static double calculateEAM(double[][] matrixOne, double[][] matrixTwo) {
        int M = matrixOne.length;
        int N = matrixOne[0].length;
        double erroAbsMed = 0;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                erroAbsMed += Math.abs(matrixOne[i][j] - matrixTwo[i][j]);
            }
        }

        return erroAbsMed / (M * N);
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

    public static boolean testCalculateEAM(double[][] A, double[][] Ak, double expectedEAM) {
        // Chama a função que calcula o EAM
        double calculatedEAM = calculateEAM(A, Ak);

        // Compara o valor calculado com o esperado
        return calculatedEAM == expectedEAM;
    }

    public static double[] averageColoumn(double[][] matrix) {
        int rows = matrix.length;
        int columns = matrix[0].length;
        double[] averageColoumn = new double[columns];
        for (int i = 0; i < columns; i++) {
            double sum = 0;
            for (int j = 0; j < rows; j++) {
                sum += matrix[j][i];
            }
            averageColoumn[i] = sum / rows;
        }
        return averageColoumn;
    }


    public static boolean testColunaMedia(double[] matrix, double[] expectedColMedia) {
        return Arrays.equals(matrix, expectedColMedia);
    }

    public static double[][] deviationMatrix(double[][] matrix, double[] averageColumn) {
        int columns = matrix[0].length;
        int rows = matrix.length;
        double[][] deviationMatrix = new double[rows][columns];

        for (int j = 0; j < columns; j++) {
            for (int i = 0; i < rows; i++) {
                //! Corrigido para usar o índice correto em averageColumn
                deviationMatrix[i][j] = matrix[i][j] - averageColumn[j];
            }
        }
        return deviationMatrix;
    }


    public static boolean testCalDesvios(double[][] matrix, double[][] expectedDesvios) {
        return Arrays.equals(matrix, expectedDesvios);
    }

    public static double[][] covariances (double[][] matrixA) {
        double[][] ATranspose = transposeMatriz(matrixA);
        double[][] AxATranspose = multiplyMatrix(matrixA,ATranspose);
        return multiplyMatrixByScalar(AxATranspose,1.0/matrixA.length);
    }

    public static boolean testCovariancia(double[][] C, double[][] expectedC) {
        return C == expectedC;
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

    public static double[][] vetoresPropriosAxAT(double[][] A,double[][] AT) {
        double[][] AxAT = multiplyMatrix(A,AT);
        EigenDecomposition eigenDecomposition = decomposeMatrix(AxAT);
        RealMatrix V = eigenDecomposition.getV();
        double[][] ui = V.getData();
        return ui;
    }

    public static double[][] valoresPropriosAxAT(double[][] A,double[][] AT) {
        double[][] AxAT = multiplyMatrix(A,AT);
        EigenDecomposition eigenDecomposition = decomposeMatrix(AxAT);
        RealMatrix D = eigenDecomposition.getD();
        double[][] dArray = D.getData();
        adjustPrecision(dArray);
        return dArray;
    }

    //ignora valores negativos muiiiito baixos, para nao obter -0.0
    private static void adjustPrecision (double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (Math.abs(matrix[i][j]) < 0.001) {
                    matrix[i][j] = 0.0;
                }
            }
        }
    }

    public static double[][] valoresPropriosC(double[][] valoresPropriosAxAT, int N) {
        double[][] lambdai = multiplyMatrixByScalar(valoresPropriosAxAT,1.0/N);
        return lambdai;
    }

    public static double[][] normalizeVectors(double[][] eigenVectors) {
        int rows = eigenVectors.length;
        int cols = eigenVectors[0].length;
        double[][] normalized = new double[rows][cols];

        for (int i = 0; i < cols; i++) {
            double norm = 0;
            for (int j = 0; j < rows; j++) {
                norm += eigenVectors[j][i] * eigenVectors[j][i];
            }
            norm = Math.sqrt(norm);

            for (int j = 0; j < rows; j++) {
                normalized[j][i] = eigenVectors[j][i] / norm;
            }
        }
        return normalized;
    }



    private static double calculateNorm(double[] vector) {
        double sum = 0;
        for (int i = 0; i < vector.length; i++) {
            sum = sum + vector[i] * vector[i];
        }
        return Math.sqrt(sum);
    }

    public static EigenDecomposition decomposeMatrix(double[][] matrixToDecompose) {
        Array2DRowRealMatrix matrix = new Array2DRowRealMatrix(matrixToDecompose);
        EigenDecomposition eigenDecomposition = new EigenDecomposition(matrix);

        return eigenDecomposition;
    }

    public static double[][] transposeMatriz(double[][] matrix) {
        double[][] matrixTransposed = new double[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                matrixTransposed[j][i] = matrix[i][j];
            }
        }
        return matrixTransposed;
    }

    public static double[][] multiplyMatrix(double[][] matrizLeft, double[][] matrizRight) {
        double[][] matrixResult = new double[matrizLeft.length][matrizRight[0].length];
        for (int i = 0; i < matrizLeft.length; i++) {
            for (int j = 0; j < matrizRight[0].length; j++) {
                for (int k = 0; k < matrizRight.length; k++) {
                    matrixResult[i][j] += matrizLeft[i][k] * matrizRight[k][j];
                }
            }
        }
        return matrixResult;
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