import org.la4j.Matrix;
import org.la4j.decomposition.EigenDecompositor;

public class Andre {
    // Aproxima uma matriz com a decomposição própria
    public static Matrix approximateMatrix(Matrix A, int k) {
        // Valida o tamanho da matriz
        validateMatrixSize(A);

        // Realiza a decomposição própria ao utilizar a lib la4j
        EigenDecompositor decompositor = new EigenDecompositor(A);
        Matrix[] eigenDecomposition = decompositor.decompose();

        // Os valores próprios estão numa matriz diagonal
        Matrix D = eigenDecomposition[0];
        // Os vetores próprios estão no espaço das colunas da segunda matriz
        Matrix P = eigenDecomposition[1];

        // Extrai os k maiores valores próprios e os seus vetores correspondentes
        Matrix Dk = D.copy();
        Matrix Pk = P.copy();

        // Zera os valores próprios e os vetores associados que não estão no top-k
        for (int i = k; i < D.columns(); i++) {
            Dk.set(i, i, 0); // Define os valores próprios restantes como zero
            for (int j = 0; j < P.rows(); j++) {
                Pk.set(j, i, 0); // Zera as colunas correspondentes nos vetores próprios
            }
        }

        // Reconstrói a matriz aproximada Ak = Pk * Dk * Pk^T
        Matrix PkDk = Pk.multiply(Dk);
        return PkDk.multiply(Pk.transpose());
    }

    // Método para validar o tamanho da matriz
    private static void validateMatrixSize(Matrix A) {
        int rows = A.rows();
        int columns = A.columns();
        if (rows > 256 || columns > 256) {
            throw new IllegalArgumentException(
                    "Erro: A matriz ultrapassa a dimensão máxima permitida de 256x256. " +
                            "Dimensão atual: " + rows + "x" + columns
            );
        }
    }
}
