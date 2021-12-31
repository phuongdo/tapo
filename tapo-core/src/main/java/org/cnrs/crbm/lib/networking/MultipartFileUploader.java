package org.cnrs.crbm.lib.networking;

import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Created by pdoviet on 2/5/2015.
 */
public class MultipartFileUploader {

    public static void main(String[] args) {
        String charset = "UTF-8";
//        File uploadFile1 = new File("e:/Test/PIC1.JPG");
//        File uploadFile2 = new File("e:/Test/PIC2.JPG");
        String requestURL = "http://protein.bio.unipd.it/raphael/Raphael.jsp";

        try {
            MultipartUtility multipart = new MultipartUtility(requestURL, charset);

            multipart.addHeaderField("User-Agent", "CodeJava");
            multipart.addHeaderField("Test-Header", "Header-Value");

            multipart.addFormField("pdbcode", "1lxa");
            multipart.addFormField("chain", "A");

//            multipart.addFilePart("fileUpload", uploadFile1);
//            multipart.addFilePart("fileUpload", uploadFile2);

            List<String> response = multipart.finish();

            // System.out.println("SERVER REPLIED:");

            for (String line : response)
                if (line.contains("URL=http://protein.bio.unipd.it/raphael/work")) {
                    System.out.println(line);
                    Document doc = Jsoup.parse("<html><head>"+line+"</head></html>");
                    String content = doc.getElementsByTag("meta").get(0).attr("content");
                    System.out.println(content.split("=")[1]);
                    break;
                }
        } catch (IOException ex) {
            System.err.println(ex);
        }
    }
}