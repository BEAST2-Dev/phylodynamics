package beast.phylodynamics.epidemiology;

import beast.core.Operator;
import beast.core.Description;
import beast.core.Input;
import beast.core.Valuable;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;


/**
 * User: Denise
 * Date: Jun 29, 2012
 * Time: 7:19:46 PM
 */
@Description("Operator that simulates new S(E)IR trajectory")
public class HybridTrajectorySingleLocationSimulator extends Operator {


    public Input<HybridSEIREpidemic> SIR =
            new Input<HybridSEIREpidemic>("SIR", "SIR trajectory calculation node - simulates dS", Input.Validate.REQUIRED);

    public Input<Valuable> S0_input =
            new Input<Valuable>("S0", "The numbers of susceptible individuals", Input.Validate.REQUIRED);

    public Input<RealParameter> dS_input =
            new Input<RealParameter>("dS", "dS vector containing the changes in numbers of susceptibles per location", Input.Validate.REQUIRED);
    public Input<RealParameter> dR_input =
            new Input<RealParameter>("dR", "dR vector containing the changes in numbers of susceptibles per location", Input.Validate.REQUIRED);

    public Input<Valuable> birth =
            new Input<Valuable>("birth", "birth rate vector with rate per location",  Input.Validate.REQUIRED);
    public Input<Valuable> death =
            new Input<Valuable>("death", "death rate vector with rate per location",  Input.Validate.REQUIRED);
    public Input<Valuable> sampling =
            new Input<Valuable>("sampling", "sampling rate vector with rate per location",  Input.Validate.REQUIRED);


    public Input<Tree> m_tree =
            new Input<Tree>("tree", "The phylogenetic tree being estimated",  Input.Validate.REQUIRED);
    public Input<RealParameter> origin =
            new Input<RealParameter>("origin", "The origin of infection (> treeheight)", Input.Validate.REQUIRED);

    int m;
    Integer S0;
    int ntaxa;

    Tree tree;
    double T;
    HybridSEIREpidemic current;

    Double b;
    Double[] d;
    Double[] s;
    Boolean birthChanges;
    Boolean deathChanges;
    Boolean samplingChanges;


    @Override
    public void initAndValidate() throws Exception {

        S0 = (int) (S0_input.get().getArrayValue());

        if (birth.get() != null && death.get() != null && sampling.get() != null){

            b = birth.get().getArrayValue();

            int dim = death.get().getDimension();
            if (dim != sampling.get().getDimension()) throw new RuntimeException("Error: Death and sampling must have equalt dimensions!");
            d = new Double[dim];
            s = new Double[dim];

            for (int i = 0; i<dim; i++){
                d[i] = death.get().getArrayValue(i);
                s[i] = sampling.get().getArrayValue(i);
            }

        }

        else{
            throw new RuntimeException("Either specify birthRate, deathRate and samplingRate OR specify R0, becomeUninfectiousRate and samplingProportion!");
        }


        m = SIR.get().Nsamples.get();

        tree = m_tree.get();
        T = origin.get().getValue();

        ntaxa = tree.getLeafNodeCount();

        // initialize trajectory
        current =  SIR.get();
        if ( !current.initTraj(b, 0., d, s, T, ntaxa, 100, m, current.times))
            throw new RuntimeException("Could not find suitable trajectory. Please try different epi parameters!");

        dS_input.get().assignFromWithoutID(new RealParameter(current.dS));
        dR_input.get().assignFromWithoutID(new RealParameter(current.dR));
    }




    @Override
    public double proposal() {

        S0 = (int) (S0_input.get().getArrayValue());


        b = birth.get().getArrayValue();

        int dim = death.get().getDimension();
        if (dim != sampling.get().getDimension()) throw new RuntimeException("Error: Death and sampling must have equalt dimensions!");
        d = new Double[dim];
        s = new Double[dim];

        for (int i = 0; i<dim; i++){
            d[i] = death.get().getArrayValue(i);
            s[i] = sampling.get().getArrayValue(i);
        }

        T =  origin.get().getValue();

        // simulate trajectory
        current =  SIR.get();
        if ( !current.refresh(S0,0, b, 0., d, s, T, ntaxa, 1, m, current.times) /* && Randomizer.nextDouble()>.75*/){
            return Double.NEGATIVE_INFINITY;
        }

        for (int i =0; i < m-1; i++){
            dS_input.get().setValue(i, current.dS[i]);
            dR_input.get().setValue(i, current.dR[i]);
        }

        return 0.;
    }

    public static void main(String[] args){
    }
}
