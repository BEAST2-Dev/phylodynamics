package beast.evolution.operators;


import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.StateNode;
import beast.core.parameter.Parameter;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;
import beast.evolution.tree.Tree;
import beast.evolution.tree.Node;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;


/**
 * User: Denise
 * Date: Jul 3, 2013
 * Time: 12:34:50 PM
 */

@Description("This element represents an operator that scales multiple parameters in the same direction. " +
        "Each operation involves selecting a scale uniformly at random between scaleFactor and 1/scaleFactor.")

public class MultiScaleOperator extends Operator {

    public Input<Double> m_scaleFactor = new Input<Double>("scaleFactor",
            "magnitude factor used for scaling", Validate.REQUIRED);
    public Input<List<StateNode>> m_parameters = new Input<List<StateNode>>("parameter",
            "zero or more items to scale upwards", new ArrayList<StateNode>());
    public Input<Boolean> m_bOptimise = new Input<Boolean>("optimise", "flag to indicate that the scale factor is automatically changed in order to acheive a good acceptance rate (default true)", true);
    public Input<Boolean> elementWise = new Input<Boolean>("elementWise", "flag to indicate that the scaling is applied to a random index in multivariate parameters (default false)", false);
    public Input<Boolean> m_pRootOnly = new Input<Boolean>("rootOnly", "scale root of a tree only, ignored if tree is not specified (default false)", false);

    double m_fScaleFactor;

    @Override
    public void initAndValidate() throws Exception {
        m_fScaleFactor = m_scaleFactor.get();
        // sanity checks
        if (m_parameters.get().size() == 0) {
            System.err.println("WARNING: At least one item must be specified");
        }

        if (elementWise.get() && m_pRootOnly.get()) throw new RuntimeException("elementwise scaling no implemented for rootonly tree scaling");
    }

    // todo: check if HRs are right! 

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal
     *         should not be accepted
     */
    @Override
    public final double proposal() {

        final double scale = (m_fScaleFactor + (Randomizer.nextDouble() * ((1.0 / m_fScaleFactor) - m_fScaleFactor)));
        int goingUp = 0;

        try{

            if (elementWise.get()) {
                int size = 0;
                for (StateNode up : m_parameters.get()) {
                    if (size == 0) size = up.getDimension();
                    if (size > 0 && up.getDimension() != size) {
                        throw new RuntimeException("elementWise=true but parameters of differing lengths!");
                    }
                    goingUp += 1;
                }


                int index = Randomizer.nextInt(size);

                for (StateNode up : m_parameters.get()) {
                    if (up instanceof RealParameter) {

                        RealParameter p = (RealParameter) up;
                        p.setValue(p.getValue(index) * scale);
                    }
                    if (outsideBounds(up)) {
                        return Double.NEGATIVE_INFINITY;
                    }
                }

            } else {

                try {
                    for (StateNode up : m_parameters.get()) {
                        up = up.getCurrentEditable(this);

                        if (up instanceof Tree && m_pRootOnly.get()) {
                            Tree tree = (Tree) up;
                            Node root = tree.getRoot();
                            double fNewHeight = root.getHeight() * scale;
                            if (fNewHeight < Math.max(root.getLeft().getHeight(), root.getRight().getHeight())) {
                                return Double.NEGATIVE_INFINITY;
                            }
                            root.setHeight(fNewHeight);
//                            return -Math.log(scale);
                            goingUp += 1;
                        } else {
                            goingUp += up.scale(scale);
                        }

                        if (outsideBounds(up)) {
                            return Double.NEGATIVE_INFINITY;
                        }
                    }

                } catch (Exception e) {
                    // scale resulted in invalid StateNode, abort proposal
                    return Double.NEGATIVE_INFINITY;
                }
            }

        } catch (Exception e) {
            // whatever went wrong, we want to abort this operation...
            return Double.NEGATIVE_INFINITY;
        }

        return (goingUp - 2) * Math.log(scale);
    }

    private boolean outsideBounds(final StateNode node) {
        if (node instanceof Parameter<?>) {
            final Parameter<?> p = (Parameter) node;
            final Double lower = (Double) p.getLower();
            final Double upper = (Double) p.getUpper();
            final Double value = (Double) p.getValue();
            if (value < lower || value > upper) {
                return true;
            }
        }
        return false;
    }

    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (m_bOptimise.get()) {
            double fDelta = calcDelta(logAlpha);
            fDelta += Math.log(1.0 / m_fScaleFactor - 1.0);
            m_fScaleFactor = 1.0 / (Math.exp(fDelta) + 1.0);
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return m_fScaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double fValue) {
        m_fScaleFactor = fValue;
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double sf = Math.pow(m_fScaleFactor, ratio);

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > 0.40) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }
} // class UpDownOperator

