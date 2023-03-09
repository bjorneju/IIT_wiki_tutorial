import pyphi
import pandas as pd
import numpy as np


def pandas_tpm(tpm, input_labels, output_labels ):
    
    past_states = list(pyphi.utils.all_states(len(past_labels)))
    future_states = list(pyphi.utils.all_states(len(future_labels)))

    columns = pd.MultiIndex.from_tuples(future_states, names=future_labels)
    index = pd.MultiIndex.from_tuples(past_states, names=past_labels)

    return pd.DataFrame(tpm, index=index, columns=columns)


def state_by_state_tpm(network,data=False):
    
    if not data:
        tpm = pyphi.convert.sbn2sbs(network.tpm)
    else:
        tpm = data

    past_labels=future_labels=list(network.node_labels) 
    
    past_states = list(pyphi.utils.all_states(len(past_labels)))
    future_states = list(pyphi.utils.all_states(len(future_labels)))

    columns = pd.MultiIndex.from_tuples(future_states, names=future_labels)
    index = pd.MultiIndex.from_tuples(past_states, names=past_labels)

    return pd.DataFrame(tpm, index=index, columns=columns)

def state_by_node_tpm(network):
    
    tpm = pyphi.convert.to_2d(network.tpm)

    past_labels=list(network.node_labels) 
    f'State by node TPM for {"".join(past_labels)}.'
    past_states = list(pyphi.utils.all_states(len(past_labels)))
    
    columns = [l+' ON' for l in network.node_labels]
    index = pd.MultiIndex.from_tuples(past_states, names=past_labels)

    return pd.DataFrame(tpm, index=index, columns=columns)

def unit_tpm(network, unit):
    
    if unit in network.node_labels:
        label = unit
        ix = network.node_indices[network.node_labels.index(label)]
        
    elif unit in network.node_indices:
        ix = unit
        label = network.node_labels[network.node_indices.index(ix)]
    else:
        print('That unit ID seems mistaken. returning an empty variable')
        return
    
    tpm = pyphi.convert.to_2d(network.tpm)[:,ix]
    past_labels=list(network.node_labels)
    past_states = list(pyphi.utils.all_states(len(past_labels)))
    
    columns = [label+' ON']
    index = pd.MultiIndex.from_tuples(past_states, names=past_labels)

    return pd.DataFrame(tpm, index=index, columns=columns)

def unit_state_by_state_tpm(network, unit):
    
    if unit in network.node_labels:
        label = unit
        ix = network.node_indices[network.node_labels.index(label)]
        
    elif unit in network.node_indices:
        ix = unit
        label = network.node_labels[network.node_indices.index(ix)]
    else:
        print('That unit ID seems mistaken. returning an empty variable')
        return
    
    tpm = pyphi.convert.to_2d(network.tpm)[:,ix]
    tpm = np.squeeze(np.dstack((tpm,1-tpm)))
    
    past_labels=list(network.node_labels)
    past_states = list(pyphi.utils.all_states(len(past_labels)))
    future_states = list(pyphi.utils.all_states(1))

    columns = pd.MultiIndex.from_tuples(future_states, names=[label])
    index = pd.MultiIndex.from_tuples(past_states, names=past_labels)

    return pd.DataFrame(tpm, index=index, columns=columns)

def possible_inputs(network,unit_label):
    
    unit_index = network.node_labels.index(unit_label)
    
    input_units = set([
        unit 
        for purview in network.potential_purviews(pyphi.direction.Direction.CAUSE,(unit_index,))
        for unit in purview 
    ])
    
    return [network.node_labels[i] for i in input_units]


def substrate_probability_from_unit_probabilities(substrate,input_state, output_state):
    return np.prod([
        p
        if unit_output==1
        else 1-p
        for unit_output, p in zip(output_state,substrate.tpm[input_state])
    ])


def unconstrained_unit_probability(substrate,unit,output_state):
    TPM = unit_tpm(substrate,unit)
    
    if output_state==1 or output_state==(1,):
        return TPM.mean()[0]
    else:
        return 1-TPM.mean()[0]
    
    
def constrained_unit_probability(substrate,unit,input_state,output_state):
    TPM = unit_state_by_state_tpm(substrate,unit)
    if type(output_state) is int:
        output_state = (output_state,)
        
    return TPM[output_state][input_state]
    
def constrained_probability(substrate,output_state,input_state):
    return np.prod([
        constrained_unit_probability(substrate,unit,input_state,unit_output_state)
        for unit,unit_output_state in zip(substrate.node_labels,output_state)
    ])
    
def unconstrained_probability(substrate,output_state):
    return np.prod([
        unconstrained_unit_probability(substrate,unit,unit_output_state)
        for unit,unit_output_state in zip(substrate.node_labels,output_state)
    ])

def all_constrained_probability(substrate):
    P = [
        [
            constrained_probability(substrate,output_state,input_state)
            for output_state in pyphi.utils.all_states(len(substrate.node_labels))
        ]
        for input_state in reversed(list(pyphi.utils.all_states(len(substrate.node_labels))))
    ]
    return state_by_state_tpm(substrate,data=P)

def all_unconstrained_probability(substrate):
    P = [
        [
            unconstrained_probability(substrate,output_state)
            for output_state in pyphi.utils.all_states(len(substrate.node_labels))
        ]
        for input_state in reversed(list(pyphi.utils.all_states(len(substrate.node_labels))))
    ]
    return state_by_state_tpm(substrate,data=P)
