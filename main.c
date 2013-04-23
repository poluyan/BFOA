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

#define S       50
#define Sr      S/2
#define ss      0.6
#define N_ed    3
#define N_re    6
#define N_ch    20
#define N_sl    4
#define p_ed    0.25
#define d_attr  0.1
#define w_attr  0.2
#define h_rep   d_attr
#define w_rep   10.0

typedef struct Cell
{
    double vect[dimension];
    double cost;
    double fitness;
    double health;
    double step_size;
} Cell;

Cell population[S];

double space[dimension][2];
double rand_vect[dimension];
double delta[dimension];

double best = INF;
int fe_count = 0;


/* functions */
void objective_function(Cell *x);
void interaction(Cell *x);

double random_number(double a, double b);
void initialize_space(double a, double b);
void initialize_population();

void tumble_step(Cell *new_cell, Cell *current_cell);
void swim_step(Cell *new_cell, Cell *current_cell);

int compare(struct Cell *left, struct Cell *right);

void chemotaxis();
void reproduction();
void elimination_dispersal();
void optimization();

int main()
{
    srand(1);

    printf("Bacterial Foraging Optimization Algorithm\n");
    printf("Dimension: %d\n", dimension);

    initialize_space(-100.0, 100.0);
    initialize_population();
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
    qsort(population, S, sizeof(Cell), (int(*)(const void*,const void*))compare);
    int i, j;
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
    for(l = 0; l < N_ed; l++)
    {
        for(k = 0; k < N_re; k++)
        {
            for(j = 0; j < N_ch; j++)
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
