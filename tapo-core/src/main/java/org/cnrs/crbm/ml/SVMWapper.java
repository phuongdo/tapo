package org.cnrs.crbm.ml;


import libsvm.svm;
import libsvm.svm_model;
import libsvm.svm_node;
import libsvm.svm_parameter;
import org.cnrs.crbm.lib.repeats.CombineScore;

import java.util.StringTokenizer;


/**
 * Created by pdoviet on 5/30/2015.
 */
public class SVMWapper {

    public static void main(String[] args) throws Exception {


        CombineScore combineScore = new CombineScore("TM-Score=0.441;Psim-Score=0.812;CE-Score=1.173;V-Score=0.515;L-Score=0.000;CA-Score=0.000;S-Score=0.441;RA=0.003");
        // System.out.println(svmScoreTRs);
        String line = "1 1:0.0 2:0.00 3:1 4:0.00 5:0.00 6:0.00 7:0.9";

        line = combineScore.toLibSVMFormat();
        System.out.println(line);
//        String line = "1 1:0.10 2:0.10 3:0.10 4:0.10 5:0.10 6:0.10 7:0.10";
        double svmScore = SVMWapper.getInstance().getProbEstimatesScore(line);
        System.out.println(svmScore);

    }


    public double getProbEstimatesScore(String line) {

        int predict_probability = 1;
        int svm_type = svm.svm_get_svm_type(model);
        int nr_class = svm.svm_get_nr_class(model);
        double[] prob_estimates = null;

        if (predict_probability == 1) {
            int[] labels = new int[nr_class];
            svm.svm_get_labels(model, labels);
            prob_estimates = new double[nr_class];

        }

        //String line = "1 1:0.977934 2:0.96988 3:0.249379 4:0.743 7:0.365023";
        StringTokenizer st = new StringTokenizer(line, " \t\n\r\f:");
        double target = atof(st.nextToken());
        int m = st.countTokens() / 2;
        svm_node[] x = new svm_node[m];
        for (int j = 0; j < m; j++) {
            x[j] = new svm_node();
            x[j].index = atoi(st.nextToken());
            x[j].value = atof(st.nextToken());
        }
        double v;
        double svmScoreTRs = 0.0;
        if (predict_probability == 1 && (svm_type == svm_parameter.C_SVC || svm_type == svm_parameter.NU_SVC)) {
            v = svm.svm_predict_probability(model, x, prob_estimates);
//            System.out.print(v + " ");
//            for (int j = 0; j < nr_class; j++)
//                System.out.print(prob_estimates[j] + " ");
//            System.out.print("\n");
            svmScoreTRs = prob_estimates[0];


        }

        return svmScoreTRs;


    }

    private static SVMWapper instance = null;

    protected SVMWapper() {
        String modeDir = "model/TRs.model";
        try {

            this.model = svm.svm_load_model(modeDir);
            if (model == null) {
                System.err.print("can't open model file " + modeDir + "\n");
                // System.exit(1);
            }
        } catch (Exception ex) {
            System.err.print("can't open model file " + modeDir + "\n");
        }
    }


    public static SVMWapper getInstance() {
        if (instance == null) {
            instance = new SVMWapper();
        }
        return instance;
    }

    svm_model model;

    private static double atof(String s) {
        return Double.valueOf(s).doubleValue();
    }

    private static int atoi(String s) {
        return Integer.parseInt(s);
    }
}
