/*
    Bacterial Foraging Optimization Algorithm proposed by Kevin M. Passino.
    Copyright (C) 2013 Sergey Poluyan <svpoluyan@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define INF DBL_MAX
#define PI acos(-1.0)

#define dimension 10

#define S       50      /* population size */
#define Sr      S/2     /* number to split */
#define ss      0.6     /* step size */
#define N_ed    3       /* number of elimination-dispersal events */
#define N_re    6       /* number of reproduction steps */
#define N_ch    20      /* number of chemotactic steps */
#define N_sl    4       /* swim length */
#define p_ed    0.25    /* eliminate probability */
#define d_attr  0.1     /* depth of the attractant */
#define w_attr  0.2     /* width of the attractant signal */
#define h_rep   d_attr  /* height of the repellant effect */
#define w_rep   10.0    /* width of the repellant */

/* bacterium */
typedef struct Cell
{
    double vect[dimension];     /* position in search space */
    double cost;                /* objective function value */
    double fitness;             /* cost value and attractant and repellent effects */
    double health;              /* the health of bacterium */
    double step_size;           /* step in the search area */
} Cell;

Cell population[S];             /* population of bacteria */

double space[dimension][2];     /* the boundaries of the search space */
double rand_vect[dimension];    /* direction of movement after a tumble */
double delta[dimension];        /* used in the normalization of the rand_vect */

double best = INF;              /* the best solution found during the search */
int fe_count = 0;               /* number of objective function evaluations */

/* functions */

/* compute objective function */
void objective_function(Cell *x);

/* compute cell-to-cell attraction and repelling effects */
void interaction(Cell *x);

/* generate random number from a to b */
double random_number(double a, double b);

/* set the bounds values for search space */
void initialize_space(double a, double b);

/* distribute the population within the search space */
void initialize_population();

/* tumble current_cell, one step in a random direction */
void tumble_step(Cell *new_cell, Cell *current_cell);

/* swim step of current_cell in a rand_vect direction */
void swim_step(Cell *new_cell, Cell *current_cell);

/* function that compares two Cell objects by health value */
int compare(struct Cell *left, struct Cell *right);

/* tumble and swim each member in the population */
void chemotaxis();

/* split the bacteria */
void reproduction();

/* elimination and dispersal event */
void elimination_dispersal();

/* run an algorithm */
void optimization();

int main()
{
    srand(1);

    printf("Bacterial Foraging Optimization Algorithm\n");
    printf("Dimension: %d\n", dimension);

    /* search space [-100, 100]^dimension */
    initialize_space(-100.0, 100.0);
    /* random initialization within the search space */
    initialize_population();
    /* minimization of objective function */
    optimization();

    return 0;
}
void objective_function(Cell *x)
{
    double rez = 0.0;
    fe_count++;

    /* Sphere Function */
    int i;
    for(i = 0; i < dimension; i++)
        rez += pow(x->vect[i], 2.0);

    x->cost = rez;

    if(x->cost < best)
        best = x->cost;
}
double random_number(double a, double b)
{
    return ((((double)rand())/((double)RAND_MAX) )*(b-a) + a);
}
void initialize_space(double a, double b)
{
    int i;
    for(i = 0; i < dimension; i++)
    {
        space[i][0] = a;
        space[i][1] = b;
    }
}
int compare(struct Cell *left, struct Cell *right)
{
    if( left->health < right->health)
        return -1;
    if (left->health > right->health)
        return 1;
    return 0;
}
void initialize_population()
{
    /* randomly distribute the initial population */
    int i, j;
    for(i = 0; i < S; i++)
    {
        for(j = 0; j < dimension; j++)
        {
            population[i].vect[j] = random_number(space[j][0], space[j][1]);
        }
        objective_function(&population[i]);
        population[i].fitness = 0.0;
        population[i].health = 0.0;
        population[i].step_size = ss;
    }
}
void elimination_dispersal()
{
    int i, j;
    for(i = 0; i < S; i++)
    {
        /* simply disperse bacterium to a random location on the search space */
        if(random_number(0.0,1.0) < p_ed)
        {
            for(j = 0; j < dimension; j++)
            {
                population[i].vect[j] = random_number(space[j][0],space[j][1]);
            }
            objective_function(&population[i]);
        }
    }
}
void reproduction()
{
    /* sort the population in order of increasing health value */
    qsort(population, S, sizeof(Cell), (int(*)(const void*,const void*))compare);
    int i, j;
    /* Sr healthiest bacteria split into two bacteria, which are placed at the same location */
    for(i = S-Sr, j = 0; j < Sr; i++, j++)
    {
        population[i] = population[j];
    }
    for(i = 0; i < S; i++)
    {
        population[i].health = 0.0;
    }
}
void interaction(Cell *x)
{
    int i, j;
    double attract = 0.0, repel = 0.0, diff = 0.0;
    for(i = 0; i < S; i++)
    {
        diff = 0.0;
        for(j = 0; j < dimension; j++)
        {
            diff += pow(x->vect[j] - population[i].vect[j], 2.0);
        }
        attract += -1.0*d_attr*exp(-1.0*w_attr*diff);
        repel += h_rep*exp(-1.0*w_rep*diff);
    }
    /* this produces the swarming effect */
    x->fitness = x->cost + attract + repel;
}
void tumble_step(Cell *new_cell, Cell *current_cell)
{
    int i;
    double a = -1.0, b = 1.0, temp1 = 0.0, temp2 = 0.0;
    for(i = 0; i < dimension; i++)
    {
        delta[i] = random_number(a, b);
        temp1 += pow(delta[i], 2.0);
    }
    temp2 = sqrt(temp1);
    for(i = 0; i < dimension; i++)
    {
        rand_vect[i] = delta[i]/temp2;
        new_cell->vect[i] = current_cell->vect[i] + current_cell->step_size*rand_vect[i];
        /* there is no need to perform search outside of the given bounds */
        if(new_cell->vect[i] < space[i][0])
            new_cell->vect[i] = space[i][0];
        if(new_cell->vect[i] > space[i][1])
            new_cell->vect[i] = space[i][1];
    }
}
void swim_step(Cell *new_cell, Cell *current_cell)
{
    int i;
    for(i = 0; i < dimension; i++)
    {
        new_cell->vect[i] = new_cell->vect[i] + current_cell->step_size*rand_vect[i];
        /* there is no need to perform search outside of the given bounds */
        if(new_cell->vect[i] < space[i][0])
            new_cell->vect[i] = space[i][0];
        if(new_cell->vect[i] > space[i][1])
            new_cell->vect[i] = space[i][1];
    }
}
void chemotaxis()
{
    double Jlast;
    Cell new_cell;
    int i, j, m;
    for(i = 0; i < S; i++)
    {
        interaction(&population[i]);
        Jlast = population[i].fitness;
        tumble_step(&new_cell, &population[i]);
        objective_function(&new_cell);
        interaction(&new_cell);
        for(j = 0; j < dimension; j++)
            population[i].vect[j] = new_cell.vect[j];
        population[i].cost = new_cell.cost;
        population[i].fitness = new_cell.fitness;
        population[i].health += population[i].fitness;
        for(m = 0; m < N_sl; m++)
        {
            if(new_cell.fitness < Jlast)
            {
                Jlast = new_cell.fitness;
                swim_step(&new_cell, &population[i]);
                objective_function(&new_cell);
                interaction(&new_cell);
                for(j = 0; j < dimension; j++)
                    population[i].vect[j] = new_cell.vect[j];
                population[i].cost = new_cell.cost;
                population[i].fitness = new_cell.fitness;
                population[i].health += population[i].fitness;
            }
            else break;
        }
    }
}
void optimization()
{
    int l, k, j;
    for(l = 0; l < N_ed; l++)           /* Elimination-dispersal loop */
    {
        for(k = 0; k < N_re; k++)       /* Reproduction loop */
        {
            for(j = 0; j < N_ch; j++)   /* Chemotaxis loop */
            {
                chemotaxis();
                printf("best=%e , fe_count=%d\n", best, fe_count);
            }
            reproduction();
        }
        elimination_dispersal();
    }
    printf("\nbest found value: %e, number of function evaluations: %d\n", best, fe_count);
}
