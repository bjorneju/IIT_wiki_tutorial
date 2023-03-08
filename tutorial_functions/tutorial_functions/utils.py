import pyphi
import pandas as pd
import numpy as np


def state_by_state_tpm(network, caption=None):
    
    tpm = pyphi.convert.sbn2sbs(network.tpm)

    past_labels=future_labels=list(network.node_labels) 
    
    past_states = list(pyphi.utils.all_states(len(past_labels)))
    future_states = list(pyphi.utils.all_states(len(future_labels)))

    columns = pd.MultiIndex.from_tuples(future_states, names=future_labels)
    index = pd.MultiIndex.from_tuples(past_states, names=past_labels)

    TPM = pd.DataFrame(tpm, index=index, columns=columns)
    (TPM.style.background_gradient(cmap="gray_r", low=0, high=1, axis=None)
        .format(precision=2)
        .set_caption(
            f'State by node TPM for {"".join(past_labels)}.'
            if caption == None
            else caption
    ))
    return TPM
    

def state_by_node_tpm(network, caption=None):
    
    tpm = pyphi.convert.to_2d(network.tpm)

    past_labels=list(network.node_labels) 
    
    past_states = list(pyphi.utils.all_states(len(past_labels)))
    
    columns = [l+' ON' for l in network.node_labels]
    index = pd.MultiIndex.from_tuples(past_states, names=past_labels)

    TPM = pd.DataFrame(tpm, index=index, columns=columns)
    (TPM.style.background_gradient(cmap="gray_r", low=0, high=1, axis=None)
        .format(precision=2)
        .set_caption(
            f'State by node TPM for {"".join(past_labels)}.'
            if caption == None
            else caption
    ))
    return TPM

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

    TPM = pd.DataFrame(tpm, index=index, columns=columns)
    (TPM.style.background_gradient(cmap="gray_r", low=0, high=1, axis=None)
        .format(precision=2)
        .set_caption(
            f'Unit TPM for {label}.'
    ))
    return TPM


def possible_inputs(network,unit_label):
    
    unit_index = network.node_labels.index(unit_label)
    
    input_units = set([
        unit 
        for purview in network.potential_purviews(pyphi.direction.Direction.CAUSE,(unit_index,))
        for unit in purview 
    ])
    
    return [network.node_labels[i] for i in input_units]