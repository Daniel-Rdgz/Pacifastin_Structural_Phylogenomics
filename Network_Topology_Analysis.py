#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Network Topology Analysis of Domain Architectures.

Constructs co-occurrence networks of Pfam domains associated with
Pacifastin cores to visualize the 'Hub' vs 'Solitary' transition.
"""

import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd

def build_network(df_domains, lineage_name):
    G = nx.Graph()
    
    # Add Central Node
    G.add_node("Pacifastin", type="Core", size=2000)
    
    # Add accessory domains
    for _, row in df_domains.iterrows():
        domain = row['Accessory_Domain']
        if domain != "Pacifastin":
            if not G.has_node(domain):
                G.add_node(domain, type="Accessory", count=0)
            
            G.nodes[domain]['count'] += 1
            
            # Add edge or increase weight
            if G.has_edge("Pacifastin", domain):
                G[ "Pacifastin"][domain]['weight'] += 1
            else:
                G.add_edge("Pacifastin", domain, weight=1)
                
    return G

# Note: Visualization code using nx.draw_networkx goes here
# Ensure to export to Gephi or Cytoscape format for high-quality rendering if needed:
# nx.write_gexf(G, f"results/{lineage_name}_network.gexf")