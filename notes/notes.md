# Notes

- Le fait de traiter la série de temps comme un stack 3D n'est pas une bonne idée.
  Cela détruit de l'information utile quand l'objet disparait sur la frame suivante et quand
  l'objet ne fait que bouger, on observe la formation de deux spots sur les frames où il a 
  été propagé.