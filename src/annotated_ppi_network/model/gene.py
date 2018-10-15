# -*- coding: utf-8 -*-

"""This module contains the class Gene."""


class Gene:
    """Encapsulates a gene and its attributes."""

    def __init__(self,
                 entrez_id: str,
                 log2_fold_change: float = 0.0,
                 symbol: str = "",
                 padj: float = 1.0):
        """Construct a Gene object.

        :param str entrez_id: Entrez id for the gene.
        :param float log2_fold_change: Log fold change of the gene.
        :param str symbol: Symbol of the gene.
        :param float padj: Adjusted p-value.
        """
        self.entrez_id = entrez_id
        self.log2_fold_change = log2_fold_change
        self.symbol = symbol
        self.padj = padj

    def __repr__(self):
        """Return a representation of the Gene object.

        :return str: A string that includes the information on HGNC symbol,
         log fold change and adjusted p-value of the gene.
        """
        return "HGNC: {symbol} L2FC: {log2_fold_change} PADJ: {padj}".format(
            symbol=self.symbol,
            log2_fold_change=self.log2_fold_change,
            padj=self.padj
        )
