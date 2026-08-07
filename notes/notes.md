# Notes

- Pour essayer de se débarasser du bruit, on commence par faire un filtre médian. 
  Dans la mesure où on est en 2D avec des objets qui sont censés peu bouger, utilisation
  d'un median 3D pour aggreger de l'info des frames d'avant et après.

| N Frames | Temps |
| -------- | ----- |
| 1        | 13.1  |
| 3        | 31.5  |
| 5        | 47.6  |

