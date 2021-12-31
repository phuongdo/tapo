package org.cnrs.crbm.ml;

import net.sf.javaml.classification.Classifier;
import net.sf.javaml.classification.KNearestNeighbors;
import net.sf.javaml.classification.bayes.NaiveBayesClassifier;
import net.sf.javaml.classification.evaluation.CrossValidation;
import net.sf.javaml.classification.evaluation.PerformanceMeasure;
import net.sf.javaml.core.Dataset;

import net.sf.javaml.core.DenseInstance;
import net.sf.javaml.core.Instance;
import net.sf.javaml.tools.data.FileHandler;
import net.sf.javaml.tools.weka.WekaClassifier;
import weka.classifiers.bayes.NaiveBayes;
import weka.classifiers.functions.SMO;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Random;
import java.util.Map;

/**
 * Created by pdoviet on 7/9/2015.
 */
public class MLDemo {

    public static void main(String[] args) throws IOException {


             /* Load data */
        //Dataset data = FileHandler.loadDataset(new File("data/iris.data"), 4, ",");
        Dataset data = FileHandler.loadDataset(new File("data/resultLargeScale/scopeDB20.train.txt"), 0, "\t");
        /* Create Weka classifier */
        //NaiveBayes smo = new NaiveBayes();
        SMO smo = new SMO();
        /* Wrap Weka classifier in bridge */
        Classifier javamlsmo = new WekaClassifier(smo);
        /* Initialize cross-validation */
        CrossValidation cv = new CrossValidation(javamlsmo);
        /* Perform cross-validation */
        Map<Object, PerformanceMeasure> pm = cv.crossValidation(data);
        /* Output results */
        System.out.println(pm);

//        Dataset data = FileHandler.loadDataset(new File("data/resultLargeScale/scopeDB20.train.txt"), 0, "\t");
//        Classifier classifier = new NaiveBayesClassifier(true, true, true);
//        classifier.buildClassifier(data);
//        double[] values = new double[]{0.761, 0, 0, 0, 0, 0};
//        Instance instance = new DenseInstance(values);
//        Object predictedClassValue = classifier.classify(instance);
//        System.out.println(predictedClassValue);
        // serialize model
//        ObjectOutputStream oos = new ObjectOutputStream(
//                new FileOutputStream("model/scopeClassifier.model"));
//        oos.writeObject(classifier);
//        oos.flush();
//        oos.close();


//        CrossValidation cv = new CrossValidation(classifier);
//        Map<Object, PerformanceMeasure> r = cv.crossValidation(data, 5, new Random(25));
//        System.out.println("Accuracy=" + r.get("a").getAccuracy());

    }
}
