from inspect import getmembers
from .library import (Equation, Observable, ConstantTerm, DerivativeOrder, latexify, multiple_replace, latex_replacements, LibraryTerm, cast_to_term, ES_safe_prod, LibraryPrime, partition)
import cadabra2
import re
import unicodedata
from cdb.utils.indices import symmetrise_slots, get_free_indices, get_dummy_indices, get_indices
from cdb.utils.node import n_indices
import itertools

#reverse lookup for latex replacements
txt_replacements = {latex_replacements[plain]:plain for plain in latex_replacements}
txt_replacements["\\cdot"] = " ·"
txt_replacements["-"] = " + -"

lowercase_greek_letters = "αβγδεζηθικλμνξοπρστυφχψω"



"""
What: Checks whether an expression is entirely symmetric or antisymmetric
Return: Boolean if the requested symmetry is true

ex: Cadabra expression
antisym: true if checking antisymmetric 
identities: Map from Cadabra expression to Cadabra expression

Change to slot_asym, needs deep = False
TBD
"""
def check_full_symmetry(ex: cadabra2.Ex, antisym: bool):
    #create copy and perform projection as Sym/AntiSym is linear
    return 0

"""
What: Cannonicalizes equation (accounting for symmetries, applies young tableau if requested is in development)
Return: Canonicalized Equation

eq: Equation
obs: List of observables (might not be needed)
str_to_ex: Map of String Representation of Primes -> Prime
exstr_to_prime: Map from Cadabra expression primes to primes
"""    
def canonicalise_and_simplify_eq(eq: Equation, obs: list, str_to_ex: dict, exstr_to_prime: dict) -> Equation:
    ex = convert_eq_to_ex(eq, obs, str_to_ex)
    cadabra2.canonicalise(ex) #canonicalise expression
    return convert_ex_to_eq(ex, exstr_to_prime)

"""
What: Turns Prime into Cadabra Expression
Return: Prime Expression

p: Library Prime 
"""
def convert_prime_to_ex(p: LibraryPrime) -> cadabra2.Ex:
    #indice pool
    cadabra2.Indices(cadabra2.Ex(r"\alpha,\beta,\gamma,\delta,\epsilon,\zeta,\eta,\theta,\iota,\kappa,\lambda,\mu,\nu,\xi,o,\pi,\rho,\sigma,\tau,\upsilon,\phi,\chi,\psi,\omega,t"), cadabra2.Ex(r'vector'))
    
    #Derivative rather than partial because of implicit combining of partial derivatives
    cadabra2.Derivative(cadabra2.Ex(r"\partial{#}"))
    cadabra2.SelfCommuting(cadabra2.Ex(r"\partial{#}"))
    prime_str = latexify(str(p))
    #Wrap exponentials
    prime_str = re.sub(r'\^(\d+)', r'^{\1}', prime_str)
    prime_list = re.split(' ', prime_str)
    prime_term = ""
    brace_count = 0
    #build prime string by adding braces
    for value in prime_list:
        if '\\partial' in value:
            prime_term += value
            prime_term += '{'
            brace_count += 1
        else:
            prime_term += value
    #finish wrapping braces
    for i in range(brace_count):
        prime_term += '}'
    return cadabra2.Ex(fr"{prime_term}")

"""
Represent [C_0, ... ,C_n] as tensor symmetric on all indices of 0...n
"""
def convert_discrete_prime_to_ex(p: LibraryPrime) -> cadabra2.Ex:
    #need way to represent coarse grained product
    pass
    
"""
What: Turns Equation into Cadabra Expression
Return: Cadabra Expression

eq: Equation
obs: list of observables
str_to_ex: Map of String Representation of Primes -> Prime
"""
def convert_eq_to_ex(eq: Equation, obs: list, str_to_ex: dict) -> cadabra2.Ex:
    #Indici field
    cadabra2.Indices(cadabra2.Ex(r"\alpha,\beta,\gamma,\delta,\epsilon,\zeta,\eta,\theta,\iota,\kappa,\lambda,\mu,\nu,\xi,o,\pi,\rho,\sigma,\tau,\upsilon,\phi,\chi,\psi,\omega,t"), cadabra2.Ex(r'vector'))
    #Derivative rather than partial because of implicit combining of partial derivatives
    cadabra2.Derivative(cadabra2.Ex(r"\partial{#}"))
    cadabra2.SelfCommuting(cadabra2.Ex(r"\partial{#}"))

    #create Expression
    ex = None
    for term, coef in zip(eq.terms, eq.coeffs):
        #Each Multiplicative Term
        #Needed since *= doesnt work with expressions
        sub_exs = []
        for prime in term.primes:
            if len(sub_exs) == 0:
                sub_exs.append(str_to_ex[str(prime)])
            else:
                sub_exs.append(str_to_ex[str(prime)] * sub_exs[len(sub_exs) - 1])
        #add up primes and account for coefficient
        if ex is None:
            ex = cadabra2.Ex(fr"{str(coef)}") * sub_exs[len(sub_exs) - 1]
        else:
            ex += cadabra2.Ex(fr"{str(coef)}") * sub_exs[len(sub_exs) - 1]
    return ex
             
"""
What: Turns Cadabra Expression into Equation
Return: Equation

ex: Cadabra Expression
exstr_to_prime: Map from cadabra expression primes to primes
"""
def convert_ex_to_eq(ex: cadabra2.Ex, exstr_to_prime: dict):
    #Add checker to see if ex is 0
    
    #coefficients
    coefs = []
    #equation terms
    eq_terms = []
    
    #term by term iteration
    for child in ex.top().terms():
        #gets coef
        coefs.append(float(child.multiplier))
        term = None
        #splits into factors, which are primes by the attachment
        for c in child.factors():
            #remove scalar
            str_rep = re.sub(r'^-?\d+\.?\d*', '', str(c))
            if term == None:
                term = exstr_to_prime[str_rep]
            else:
                term = ES_safe_prod(term, exstr_to_prime[str_rep])
        term = cast_to_term(term) #turn from type LibraryPrime to type LibraryTerm
        eq_terms.append(term)
        
    return Equation(eq_terms, coefs)
        