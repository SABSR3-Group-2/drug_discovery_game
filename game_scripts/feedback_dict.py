"""
Script containing the dictionary with feedback for the user
based on their molecule choices
"""
feedback_pic50 = {
    'low': "Your molecule's pIC50 is lower than that of the target molecule. A lower binding affinity means that your drug will not bind to its target as strongly, and thus a higher dose will need to be taken for it to be effective.\n",
    'good': "Your molecule's pIC50 is acceptably close to the target value. Your molecule will bind strongly to the target, so a lower dose is needed, which can minimise side effects.\n"
}
feedback_lipophilicity = {
    'low': "Your drug's lipophilicity is too low, thus it is likely to have poor ADMET properties. Adding alkyl chains or non-polar functional groups can increase a drug's lipophilicity.\n",
    'good': "Your drug's lipophilicity is good, thus it is likely to have good ADMET properties with minimal off-target effects.\n",
    'high': "Your drug is too lipophilic. This can lead to a large number of off-target effects, increasing side effects and potential toxicity. Adding polar functional groups to your drug can help to reduce this.\n"
}

feedback_clearance = {
    'high': 'Your drug clearance is too high, which means it has a low metabolic stability. Patients will have to take large doses regularly in order to maintain a therapeutic concentration within their bloodstream, which is unpopular.\n',
    'low': "Your drug's metabolic stability is good, so it is not cleared from the body too quickly, thus very frequent dosing is not necessary.\n"
}