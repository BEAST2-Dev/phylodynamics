package beast.phylodynamics.util;
import beast.evolution.tree.Tree;

import java.awt.*;
import java.io.*;
import java.util.Random;
import java.util.Arrays;



/**
 * User: dkuh004
 * Date: Nov 1, 2010
 */
public class Stuff {



    /**
     * @param t the time in question
     * @return the index of the given time in the times array, or if the time is not in the array the index of the time
     *         next smallest
     */
    public static int index(double t, double[] times) {

        int epoch = Arrays.binarySearch(times, t);

        if (epoch < 0) {
            epoch = -epoch - 1;
        }

        return epoch;
    }


    public static double round(double d, int round){
        double factor = Math.pow(10, round);
        return Math.round(d * factor) / factor;
    }



    private Integer[] toIntegerArray(int[] array){
        int l =array.length;
        Integer[] IntegerArray = new Integer[l];
        for (int i = 0; i<l; i++){
               IntegerArray[i]=(Integer)array[i];
        }
        return IntegerArray;
    }

    public static Color[] colorArray(int n){
        Color[] color = new Color[n];
        for (int i=0; i < n; i++){
                 color[i]=Color.black;
        }
        // for more than 1 pop (and here up to 3, color the locations differently)
        if (n > 1){
            color[0]=Color.blue;
            color[1]=Color.red;
            if (n>2) color[2]=Color.orange;
        }

        return color;
    }

    public static StringBuilder addTikzTexFrame(StringBuilder builder, String xscale, String yscale){
        /*just tikz:      */
        builder.insert(0, "\\documentclass{article} \n\\usepackage{tikz} \n\\begin{document} \n \\begin{tikzpicture}[xscale = " + xscale + ", yscale = " + yscale + "] \n");
        builder.append("\n \\end{tikzpicture}\n \\end{document}");

//        /* pgfplots:      */
//        builder.insert(0, "\\documentclass{article} \n\\usepackage{tikz} \n\\usepackage{pgfplots} \n\\begin{document} \n \\begin{tikzpicture} \n \\begin{axis}[xlabel=Days,ylabel=Value] \n");
//        builder.append("\n \\end{axis} \n \\end{tikzpicture}\n \\end{document}");
        return builder;
    }


    
    public static void writeFile(String fileName, String input, Boolean open){
        try{
          PrintWriter writer = new PrintWriter(new FileWriter(fileName));
          writer.write(input);
          writer.flush();
          writer.close();
          if (open) Runtime.getRuntime().exec("open " + fileName);
        }catch(Exception e){ System.out.println(e);}

    }


     private static void treeToPDF(String filename){
          System.out.println(filename);
          String newfile = filename + ".pdf";
          String cmdline = "java -cp /Users/dkuh004/Documents/IdeaProjects/FigTree/dist/figtree-pdf.jar figtree.application.FigTreePDF" + " " + filename + " " + newfile;
          try{
          Runtime.getRuntime().exec(cmdline);
          }catch(Exception e){System.out.println(e);}
     }



    private static String getFigTreeSettings(){
           return "begin figtree;\n" +
                   "\tset appearance.backgroundColorAttribute=\"User Selection\";\n" +
                   "\tset appearance.backgroundColour=#-1;\n" +
                   "\tset appearance.branchColorAttribute=\"state\";\n" +
                   "\tset appearance.branchLineWidth=1.0;\n" +
                   "\tset appearance.foregroundColour=#-16777216;\n" +
                   "\tset appearance.selectionColour=#-2144520576;\n" +
                   "\tset branchLabels.colorAttribute=\"User Selection\";\n" +
                   "\tset branchLabels.displayAttribute=\"Branch times\";\n" +
                   "\tset branchLabels.fontName=\"Abadi MT Condensed Extra Bold\";\n" +
                   "\tset branchLabels.fontSize=8;\n" +
                   "\tset branchLabels.fontStyle=0;\n" +
                   "\tset branchLabels.isShown=false;\n" +
                   "\tset branchLabels.significantDigits=4;\n" +
                   "\tset layout.expansion=0;\n" +
                   "\tset layout.layoutType=\"RECTILINEAR\";\n" +
                   "\tset layout.zoom=0;\n" +
                   "\tset nodeBars.barWidth=4.0;\n" +
                   "\tset nodeLabels.colorAttribute=\"User Selection\";\n" +
                   "\tset nodeLabels.displayAttribute=\"Node ages\";\n" +
                   "\tset nodeLabels.fontName=\"Abadi MT Condensed Light\";\n" +
                   "\tset nodeLabels.fontSize=12;\n" +
                   "\tset nodeLabels.fontStyle=0;\n" +
                   "\tset nodeLabels.isShown=false;\n" +
                   "\tset nodeLabels.significantDigits=4;\n" +
                   "\tset polarLayout.alignTipLabels=false;\n" +
                   "\tset polarLayout.angularRange=0;\n" +
                   "\tset polarLayout.rootAngle=0;\n" +
                   "\tset polarLayout.rootLength=100;\n" +
                   "\tset polarLayout.showRoot=true;\n" +
                   "\tset radialLayout.spread=0.0;\n" +
                   "\tset rectilinearLayout.alignTipLabels=true;\n" +
                   "\tset rectilinearLayout.curvature=0;\n" +
                   "\tset rectilinearLayout.rootLength=100;\n" +
                   "\tset scale.offsetAge=0.0;\n" +
                   "\tset scale.rootAge=1.0;\n" +
                   "\tset scale.scaleFactor=1.0;\n" +
                   "\tset scale.scaleRoot=false;\n" +
                   "\tset scaleAxis.automaticScale=true;\n" +
                   "\tset scaleAxis.fontSize=8.0;\n" +
                   "\tset scaleAxis.isShown=true;\n" +
                   "\tset scaleAxis.lineWidth=1.0;\n" +
                   "\tset scaleAxis.majorTicks=1.0;\n" +
                   "\tset scaleAxis.origin=0.0;\n" +
                   "\tset scaleAxis.reverseAxis=false;\n" +
                   "\tset scaleAxis.showGrid=true;\n" +
                   "\tset scaleAxis.significantDigits=4;\n" +
                   "\tset scaleBar.automaticScale=true;\n" +
                   "\tset scaleBar.fontSize=12.0;\n" +
                   "\tset scaleBar.isShown=false;\n" +
                   "\tset scaleBar.lineWidth=1.0;\n" +
                   "\tset scaleBar.scaleRange=0.0;\n" +
                   "\tset scaleBar.significantDigits=4;\n" +
                   "\tset tipLabels.colorAttribute=\"User Selection\";\n" +
                   "\tset tipLabels.displayAttribute=\"Names\";\n" +
                   "\tset tipLabels.fontName=\"Abadi MT Condensed Light\";\n" +
                   "\tset tipLabels.fontSize=12;\n" +
                   "\tset tipLabels.fontStyle=0;\n" +
                   "\tset tipLabels.isShown=true;\n" +
                   "\tset tipLabels.significantDigits=4;\n" +
                   "\tset trees.order=false;\n" +
                   "\tset trees.orderType=\"increasing\";\n" +
                   "\tset trees.rooting=false;\n" +
                   "\tset trees.rootingType=\"User Selection\";\n" +
                   "\tset trees.transform=false;\n" +
                   "\tset trees.transformType=\"cladogram\";\n" +
                   "end;";
       }


//    public static void main(String[] args){

        //figtree.application.figTreeGraphic.
//        String[] inFiles= new String[] {"sampleTree0.tree", "wholeTree0.tree", "dyn.pdf"};
//        String out = "all.pdf";
//        joinTreePDF(inFiles, out, 1000, 1000);

        //jCompsToEPS(new JComponent[]{tree}, "test.eps");
//        String filename = "sampleTree0.tree";
//        ExtendedTreeViewer treeViewer = figtree.application.figTreeGraphic.createGraphic(100,100, filename);
//        int width = treeViewer.getContentPane().getWidth();
//        int height = treeViewer.getContentPane().getHeight();
//        Document document = new Document();
//        document.setPageSize(new com.itextpdf.text.Rectangle(width, height));
//
//        try {
//            FileOutputStream stream = new FileOutputStream(filename);
//
//            PdfWriter writer = PdfWriter.getInstance(document, stream);
//            document.open();
//            PdfContentByte cb = writer.getDirectContent();
//            PdfTemplate tp = cb.createTemplate(width, height);
//            Graphics2D g2 = tp.createGraphics(width, height);
//            tp.setWidth(width);
//            tp.setHeight(height);
//            treeViewer.getContentPane().print(g2);
//            g2.dispose();
//            tp.sanityCheck(); // all the g2 content is written to tp, not cb
//            cb.addTemplate(tp, 0, 0);
//            cb.sanityCheck();
//        } catch (Exception e) {
//            System.err.println(e.getMessage());
//        }
//        document.close();

//        try
//        {
//            FileOutputStream finalImage = new FileOutputStream(filename);
//            ExtendedTreeViewer treeviewer = figtree.application.figTreeGraphic.createGraphic(100,100, filename);
//           // Dimension d = new Dimension(treeviewer.getContentPane().getHeight(),treeviewer.getContentPane().getWidth());//tree.getSize();
//            EpsGraphics2D g = new EpsGraphics2D("", finalImage, 0, 0, 500, 1000);
//            System.out.println("g created. printing..");
//            treeviewer.getContentPane().print(g);
//            System.out.println("printed");
//
//            g.flush();
//            g.close();
//
//            finalImage.close();
//            System.out.println("File written: " + filename);
//        }
//        catch(Exception e) {System.out.println(e);}
//    }

}
