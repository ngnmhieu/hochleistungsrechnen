#include <stdio.h>
#define MAP_SIZE 3

// Definieren Sie ein enum cardd
typedef enum {
  N = 0b0001,
  S = 0b0010,
  E = 0b0100,
  W = 0b1000
  // D.h. NW=1001, NE=0101, SE=0110, SW=1010
} cardd;

// Definieren Sie ein 3x3-Array namens map, das Werte vom Typ cardd enthält
cardd map[MAP_SIZE][MAP_SIZE];

// Die Funktion set_dir soll an Position x, y den Wert dir in das Array map eintragen
// Überprüfen Sie x und y um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit
void set_dir (int x, int y, cardd dir)
{
  // Überprüft die Gültigkeit von x und y
  if (x >= MAP_SIZE || y >= MAP_SIZE)
    return;

  // Himmelsrichtung setzen wenn dir gültig ist
  switch((int)dir) {
    case N:
    case E:
    case S:
    case W:
    case N|W:
    case N|E:
    case S|E:
    case S|W:
      map[x][y] = dir;
  }
}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map (void)
{
  for (int i = 0; i < MAP_SIZE; i++) {
    for (int j = 0; j < MAP_SIZE; j++) {

      int dir = (int)map[i][j];
      switch(dir) {
        case N: printf("N"); break;
        case E: printf("E"); break;
        case S: printf("S"); break;
        case W: printf("W"); break;
        case N|W: printf("NW"); break;
        case N|E: printf("NE"); break;
        case S|E: printf("SE"); break;
        case S|W: printf("SW"); break;
        default:
          printf("0");
      }
      
      // Leerzeichen ausgeben
      if (dir == (N|W) || dir == (S|W)) {
        printf("  ");
      } else {
        printf("   ");
      }
    }
    printf("\n");
  }
}

int main (void)
{
	// In dieser Funktion darf nichts verändert werden!
	set_dir(0, 1, N);
	set_dir(1, 0, W);
	set_dir(1, 4, W);
	set_dir(1, 2, E);
	set_dir(2, 1, S);

	show_map();

	set_dir(0, 0, N|W);
	set_dir(0, 2, N|E);
	set_dir(0, 2, N|S);
	set_dir(2, 0, S|W);
	set_dir(2, 2, S|E);
	set_dir(2, 2, E|W);

	show_map();

	return 0;
}
