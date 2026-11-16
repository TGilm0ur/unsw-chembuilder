"""
Syllabus explainers shown by the 'Connect to the Syllabus' button.

Written against the NSW Chemistry Stage 6 syllabus (Year 11 content only).
Edit the text here - it is passed straight to the viewer. Light HTML (<b>, <i>, <sub>) is allowed.
"""

SYLLABUS = {
    "bonds": {
        "title": "Bond Lengths",
        "module": "Module 1: Properties and Structure of Matter  ›  Bonding",
        "question": "What binds atoms together in elements and compounds?",
        "paras": [
            "A <b>covalent bond</b> forms when two atoms share electrons. The bond length is the distance between the two nuclei where the attraction (nuclei pulling on the shared electrons) balances the repulsion (nucleus vs nucleus and electron vs electron).",
            "<b>More shared pairs, shorter bond.</b> A double bond has two shared pairs pulling the nuclei together, and a triple bond has three, so C≡C &lt; C=C &lt; C–C. Build ethane, ethene and ethyne and compare.",
            "<b>Atom size matters too.</b> Hydrogen is tiny, so X–H bonds are the shortest single bonds. Bonds to bigger atoms are longer.",
            "Shorter bonds are usually stronger bonds (more energy to break). Check the <i>Thermodynamics</i> view to see bond enthalpies for the same bonds.",
            "<i>Note:</i> lengths here come from a simple molecular-mechanics model, so treat them as good to a few per cent."
        ],
    },
    "charges": {
        "title": "Partial Charges",
        "module": "Module 1: Properties and Structure of Matter  ›  Bonding",
        "question": "Electronegativity and the spectrum of bonds from non-polar covalent to ionic.",
        "paras": [
            "<b>Electronegativity</b> is how strongly an atom attracts the shared electrons in a bond. Across this app: O &gt; N &gt; C &asymp; H (O is the most electronegative, H and C are similar).",
            "When two different atoms share electrons, the more electronegative atom pulls the pair closer and gains a small negative charge (<b>δ−</b>), leaving the other atom with a small positive charge (<b>δ+</b>). This is a <b>polar covalent bond</b>.",
            "The bigger the electronegativity difference, the bigger the partial charges. A big enough difference ends up at an ionic bond, so bonds sit on a <i>spectrum</i> from non-polar covalent to polar covalent to ionic.",
            "Try: compare the C–H bonds in methane with the O–H bonds in water. Which atoms carry the biggest charges, and why?",
            "<i>Note:</i> charges are estimated from electronegativity-based increments. They are not measured values."
        ],
    },
    "dipole": {
        "title": "Dipole Moment",
        "module": "Module 1: Properties and Structure of Matter  ›  Bonding",
        "question": "Molecular shape, polarity and the forces between molecules.",
        "paras": [
            "A <b>polar bond</b> has a small positive end and a small negative end. A whole molecule is polar only if those bond dipoles <b>do not cancel</b>. The <b>dipole moment</b> (in debye, D) adds up every bond dipole as a vector, so <b>shape</b> decides the answer.",
            "Carbon dioxide has two very polar C=O bonds, but it is linear, so they point in opposite directions and cancel (dipole ≈ 0). Water has polar O–H bonds in a bent shape, so they add up to a large dipole.",
            "<b>Polarity changes physical properties.</b> Polar molecules attract each other through dipole–dipole forces, and molecules with O–H or N–H bonds can form hydrogen bonds. Stronger intermolecular forces mean higher melting and boiling points.",
            "Try: build methane, methanol (CH<sub>3</sub>OH) and ammonia. Which are polar? Which would you expect to boil at the highest temperature?",
            "<i>Arrow convention:</i> the arrow points toward the negative end (δ−), with a + at the tail. Some textbooks and programs draw it the other way."
        ],
    },
    "ir": {
        "title": "IR Spectrum",
        "module": "Extension: not part of the Year 11 syllabus",
        "question": "How can we 'see' bonds?",
        "paras": [
            "This view goes beyond the Year 11 syllabus, but it builds on ideas you do meet: atomic mass (Module 1) and bond strength (Module 4).",
            "Bonds are not rigid sticks. They stretch and bounce like <b>springs</b>. A bond absorbs infrared light that matches its natural vibration frequency, and that shows up as a dip in the spectrum.",
            "<b>Lighter atoms and stronger bonds vibrate faster</b> (higher wavenumber). That is why O–H and C–H absorb at high wavenumbers, while heavier C–O and C–C bonds sit lower, and triple bonds sit higher than double bonds.",
            "Chemists use these 'fingerprints' to identify functional groups. Click a peak to watch that bond stretch.",
            "<i>Note:</i> this simplified model only includes stretching vibrations, so real spectra look busier."
        ],
    },
    "thermo": {
        "title": "Thermodynamics",
        "module": "Module 4: Drivers of Reactions  \u203a  Enthalpy and Hess's Law; Entropy and Gibbs Free Energy",
        "question": "How much energy does it take to break bonds, and how can enthalpy and entropy explain whether a reaction is spontaneous?",
        "paras": [
            "<b>Enthalpy:</b> breaking bonds takes in energy (endothermic) and forming bonds releases energy (exothermic). Using <b>Hess's Law</b>, the enthalpy of formation here is worked out by atomising the elements and then making your molecule's bonds. The number on each bond in the 3D view is its average bond enthalpy.",
            "<b>Entropy (S)</b> measures how spread out the energy and particles are (disorder). It rises with more gas molecules, bigger and heavier molecules, and higher temperature. Look at whether forming your molecule from its elements gives more or fewer moles of gas.",
            "<b>Gibbs free energy:</b> \u0394G = \u0394H \u2212 T\u0394S. It describes a <i>reaction</i>, so here it is worked out for the formation of your molecule from its elements (\u0394G<sub>f</sub>). A reaction is spontaneous when \u0394G is negative. Use the temperature slider to see how temperature changes the balance between \u0394H and T\u0394S.",
            "The molar mass is from the Mole Concept in Module 2 (n = m/MM). Estimates use average bond enthalpies, so expect them to differ from data-book values."
        ],
    },
}
