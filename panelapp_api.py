# Genomics England PanelApp interface

import logging
import requests
from collections import OrderedDict

__all__ = ['PanelAppQuery', 'PanelAppApi']


class PanelAppQuery(object):
    """
    Generic PanelApp query container
    """

    class FormatException(Exception):
        pass


    def __init__(self, pid, pname, pconf):
        if type(pid) != int:
            raise FormatException('PanelApp identifier: integer expected')
        if type(pname) != str:
            raise FormatException('PanelApp name: character string expected')
        if pconf not in [1, 2, 3]:
            raise FormatException('PanelApp confidence level: integer expected')
        self._pid = pid
        self._pname = pname
        self._pconf = pconf


    @property
    def id(self):
        return self._pid

    @property
    def name(self):
        return self._pname

    @property
    def confidence_level(self):
        return self._pconf


class PanelAppApi(object):
    """
    Queries Genomics England PanelApp for genes belonging in specific gene panels.

    See https://panelapp.genomicsengland.co.uk for more information
    """

    panelappBaseUrl ="https://panelapp.genomicsengland.co.uk/api/v1"

    def __init__(self, confidence_threshold=2):
        self._confid = confidence_threshold


    @staticmethod
    def get_all_panels():
        """Return dict of all GEN panels from the /panels endpoint."""
        panels = {}
        requestString = panelappBaseUrl + "/panels"
        while requestString:
            r = requests.get(requestString)
            response = r.json()
            for panel in response["results"]:
                assert panel["name"] not in panels
                panels[panel["name"]] = panel
            requestString = response["next"]
        return panels


    @staticmethod
    def get_data_for_gene_panel(panelID):
        """Get panel data (JSON) for a given panelID from /panels/{id}."""
        requestString = "/".join((panelappBaseUrl, "panels", str(panelID)))
        print("Requesting GEN {}".format(panelID))
        r = requests.get(requestString)
        return r.json()


    @staticmethod
    def gene_symbols_for_gene_panel(panelData):
        """Return dict of {hgnc, set(aliases)} from a panel."""
        gene_symbols = {}
        for gene in panelData["genes"]:
            try:
                if int(gene["confidence_level"]) < self._confid:
                    continue
            except ValueError:
                print("Error with confidence_level: {}".format(gene["confidence_level"]))
                continue
            hgnc = gene["gene_data"]["hgnc_symbol"]
            alias = gene["gene_data"]["alias"]
            gene_symbols[hgnc] = set(alias)
        return gene_symbols


    @staticmethod
    def gene_panels_and_genes(panels):
        """Return dict {panel_name:{hgnc:set(aliases)}}."""
        panel_name_gene_symbols = {}
        for panelName, panel in panels.items():
            gene_symbols_and_aliases = gene_symbols_for_gene_panel(panel)
            panel_name_gene_symbols[panelName] = _sort_dict(gene_symbols_and_aliases)
        sorted_dict = _sort_dict(panel_name_gene_symbols)
        return sorted_dict


    @staticmethod
    def _sort_dict(d):
        # Note: as of python 3.7 OrderedDict is not necessary to preserve insertion order
        sorted_dict = OrderedDict()
        for key in sorted(d):
            sorted_dict[key] = d[key]
        return sorted_dict

