package beast.phylodynamics.epidemiology;

/**
 * @author Alexei Drummond
 */

import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.RandomTree;
import beast.evolution.tree.coalescent.PopulationFunction;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * @author remco@cs.waikato.ac.nz
 */
@Description("Performs random sequence generation for a given site model. " +
        "Sequences for the leave nodes in the tree are returned as an alignment.")
public class CoalescentSimulator extends beast.core.Runnable {
    public Input<PopulationFunction> populationFunctionInput
            = new Input<PopulationFunction>("populationFunction", "the population function", Input.Validate.REQUIRED);
    public Input<Integer> ntaxaInput = new Input<Integer>("ntaxa", "the number of taxa", Input.Validate.REQUIRED);
    public Input<Integer> replicatesInput = new Input<Integer>("replicates", "the number of replicates", Input.Validate.REQUIRED);
    public Input<String> outputFileNameInput = new Input<String>(
            "outputFileName",
            "If provided, simulated trees are written to this file rather "
                    + "than to standard out.");

    /**
     * nr of categories in site model *
     */
    int ntaxa;
    /**
     * nr of states in site model *
     */
    int replicates;

    PopulationFunction populationFunction;

    /**
     * name of output file *
     */
    String outputFileName;

    TaxonSet taxa = new TaxonSet();

    @Override
    public void initAndValidate() throws Exception {

        replicates = replicatesInput.get();
        ntaxa = ntaxaInput.get();
        populationFunction = populationFunctionInput.get();
        outputFileName = outputFileNameInput.get();

        List<Taxon> taxaList = new ArrayList<Taxon>();
        for (int i = 1; i <= ntaxa; i++) {
            taxaList.add(new Taxon(i + ""));
        }
        taxa.initByName("taxon", taxaList);
    }

    @Override
    public void run() throws Exception {

        // Write output to stdout or file
        PrintStream pstream;
        if (outputFileName == null)
            pstream = System.out;
        else
            pstream = new PrintStream(outputFileName);


        for (int i = 0; i < replicates; i++) {
            RandomTree tree = new RandomTree();
            tree.initByName("taxonset", taxa, "populationModel", populationFunction);

            pstream.print(tree.toString() + "\n");
        }
        if (outputFileName != null) {
            pstream.flush();
            pstream.close();
        }
    }
} // class SequenceSimulator

