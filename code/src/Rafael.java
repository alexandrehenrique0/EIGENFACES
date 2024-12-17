public class Rafael {
    public static void main(String[] args){
        double[][] matrizVetores = {
            {1, 8, 1, 10, 9},
            {2, 3, 2, 5, 6},
            {3, 4, 3, 6, 7},
            {4, 5, 4, 7, 8},
            {8, 6, 5, 8, 9}
        };
        double[] vetorPrincipal = {1, 2, 3, 4, 5};
        double[] resultado = calculate_Euclidian_Distance(vetorPrincipal, matrizVetores);
        int min_Pos = check_Closer_Vetor(resultado);
        System.out.println("O vetor mais próximo é o vetor " + min_Pos);
    }
    public static double[] calculate_Euclidian_Distance(double[] vetorPrincipal, double[][] matrizVetores) {
        double[] resultado = new double[matrizVetores[0].length];
        for (int i = 0; i < matrizVetores[0].length; i++) {
            double soma = 0;
            for (int j = 0; j < vetorPrincipal.length; j++) {
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
}